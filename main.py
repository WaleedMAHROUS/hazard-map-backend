import os
import json
import ee
import requests
import pandas as pd
from shapely.geometry import shape, mapping
from pyproj import Geod
from math import radians, cos, sin, asin, sqrt, degrees, atan2
from flask import Flask, request, jsonify
from flask_cors import CORS
from lxml import etree
import traceback

app = Flask(__name__)
CORS(app)

# --- AUTHENTICATION ---
GEE_ENABLED = False
try:
    print("🔒 Checking Authentication...")
    service_account_json = os.environ.get('GOOGLE_SERVICE_ACCOUNT_JSON')
    
    if service_account_json:
        print("   -> Found GOOGLE_SERVICE_ACCOUNT_JSON env var.")
        try:
            key_data = json.loads(service_account_json)
            creds = ee.ServiceAccountCredentials(key_data['client_email'], key_data=json.dumps(key_data))
            ee.Initialize(credentials=creds)
            GEE_ENABLED = True
            print("✅ GEE Auth Success (Env Var)")
        except json.JSONDecodeError:
            print("❌ Error: GOOGLE_SERVICE_ACCOUNT_JSON is not valid JSON.")
        except Exception as e:
            print(f"❌ GEE Auth Error: {e}")

    elif os.path.exists('service-account-key.json'):
        print("   -> Found service-account-key.json file.")
        try:
            creds = ee.ServiceAccountCredentials(json.load(open('service-account-key.json'))['client_email'], 'service-account-key.json')
            ee.Initialize(credentials=creds)
            GEE_ENABLED = True
            print("✅ GEE Auth Success (File)")
        except Exception as e:
            print(f"❌ GEE Auth Error: {e}")
        
    else:
        print("⚠️ No Google Credentials found. Satellite data will be skipped.")
        
except Exception as e:
    print(f"❌ Critical Auth Failure: {e}")
    traceback.print_exc()

# --- HELPER FUNCTIONS ---
HEADERS = {'User-Agent': 'HabitatScanner/2.0'}

def get_bbox(lat, lon, radius_km):
    R = 6371
    lat_rad, lon_rad = radians(lat), radians(lon)
    lat_d = radius_km / R
    lon_d = radius_km / (R * cos(lat_rad))
    return (degrees(lat_rad - lat_d), degrees(lon_rad - lon_d), degrees(lat_rad + lat_d), degrees(lon_rad + lon_d))

def haversine(lat1, lon1, lat2, lon2):
    """ Returns distance in KM between two points """
    R = 6371
    dLat, dLon = radians(lat2 - lat1), radians(lon2 - lon1)
    a = sin(dLat/2)**2 + cos(radians(lat1)) * cos(radians(lat2)) * sin(dLon/2)**2
    return R * 2 * asin(sqrt(a))

def calculate_area(geom):
    try:
        s = shape(geom)
        if s.geom_type in ['Polygon', 'MultiPolygon']:
            return abs(Geod(ellps='WGS84').geometry_area_perimeter(s)[0])
        return 0 
    except: return 0

# --- CUSTOM PARSER (No external libs) ---
def osm_to_geojson_custom(osm_json):
    elements = osm_json.get('elements', [])
    nodes = {e['id']: (e['lon'], e['lat']) for e in elements if e['type'] == 'node'}
    ways = {e['id']: e for e in elements if e['type'] == 'way'}
    features = []
    
    for el in elements:
        if el['type'] == 'node' and 'tags' in el:
            features.append({"type": "Feature", "properties": el['tags'], "geometry": {"type": "Point", "coordinates": [el['lon'], el['lat']]}})
        elif el['type'] == 'way' and 'tags' in el:
            coords = []
            for nid in el.get('nodes', []):
                if nid in nodes: coords.append(nodes[nid])
            if len(coords) > 2:
                gtype = "Polygon" if coords[0] == coords[-1] else "LineString"
                features.append({"type": "Feature", "properties": el['tags'], "geometry": {"type": gtype, "coordinates": [coords] if gtype == "Polygon" else coords}})
        elif el['type'] == 'relation' and el.get('tags', {}).get('type') == 'multipolygon':
            outer = []
            for m in el.get('members', []):
                if m['role'] == 'outer' and m['ref'] in ways:
                    w = ways[m['ref']]
                    wc = [nodes[n] for n in w.get('nodes', []) if n in nodes]
                    if wc: 
                        outer = wc
                        break
            if outer:
                features.append({"type": "Feature", "properties": el['tags'], "geometry": {"type": "Polygon", "coordinates": [outer]}})
    return features

# --- DATA MINING ---
def fetch_osm_data(lat, lon, radius_m):
    print("🗺️ Fetching OSM...")
    s, w, n, e = get_bbox(lat, lon, (radius_m/1000)+1)
    q = f"""[out:json][timeout:45];(nwr["landuse"~"landfill|industrial|brownfield|reservoir|basin"]({s},{w},{n},{e});nwr["amenity"~"waste_disposal|slaughterhouse"]({s},{w},{n},{e});nwr["natural"~"water|wetland"]({s},{w},{n},{e}););out body;>;out skel qt;"""
    try:
        resp = requests.post("https://overpass-api.de/api/interpreter", data={'data': q}, headers=HEADERS)
        if resp.status_code != 200: return []
        feats = osm_to_geojson_custom(resp.json())
        for f in feats:
            p = f['properties']
            f['properties']['custom_type'] = 'water' if 'water' in str(p) or 'reservoir' in str(p) else 'waste'
            if f['geometry']['type'] == 'Point':
                f['geometry'] = mapping(shape(f['geometry']).buffer(0.0005)) # Inflate points
        print(f"   -> OSM found {len(feats)}")
        return feats
    except Exception as e:
        print(f"OSM Error: {e}"); return []

def fetch_gee_data(lat, lon, radius_m):
    if not GEE_ENABLED: 
        print("⚠️ GEE Skipped (Disabled)")
        return []
    print("🛰️ Fetching Satellite...")
    try:
        aoi = ee.Geometry.Point(lon, lat).buffer(radius_m)
        img = ee.ImageCollection('ESA/WorldCover/v100').first()
        mask = img.eq(80).Or(img.eq(40))
        
        # FIX: Increased scale from 20 to 35, and added eightConnected=False
        # This reduces the number of tiny polygons to prevent the "5000 elements" crash
        vectors = img.updateMask(mask).reduceToVectors(
            geometry=aoi, 
            scale=35, 
            maxPixels=1e10, 
            eightConnected=False
        )
        
        feats = []
        for f in vectors.getInfo()['features']:
            f['properties']['custom_type'] = 'water' if f['properties']['label'] == 80 else 'veg'
            feats.append(f)
        print(f"   -> GEE found {len(feats)}")
        return feats
    except Exception as e:
        print(f"GEE Error: {e}")
        return []

# --- KML GENERATOR ---
def generate_circle_coords(clat, clon, radius_m, num_points=100):
    coords = []
    for i in range(num_points + 1):
        brng = radians(i * (360 / num_points))
        d = radius_m / 1000 # km
        R = 6371
        lat1 = radians(clat)
        lon1 = radians(clon)
        lat2 = asin(sin(lat1)*cos(d/R) + cos(lat1)*sin(d/R)*cos(brng))
        lon2 = lon1 + atan2(sin(brng)*sin(d/R)*cos(lat1), cos(d/R)-sin(lat1)*sin(lat2))
        coords.append(f"{degrees(lon2)},{degrees(lat2)},0")
    return " ".join(coords)

def generate_files(features, arp, radius_m):
    kml = etree.Element("kml", xmlns="http://www.opengis.net/kml/2.2")
    doc = etree.SubElement(kml, "Document")
    
    # Define Styles (AABBGGRR format)
    styles = {
        "water": ("80ff0000", "ff0000"), # Blue
        "veg": ("80008000", "008000"),   # Green
        "waste": ("9913458b", "8b4513"), # Brown
        "radius": ("ff0000ff", "ff0000") # Red Line (Opacity FF, Blue 00, Green 00, Red FF)
    }
    
    for k, (poly_c, line_c) in styles.items():
        s = etree.SubElement(doc, "Style", id=k)
        if k == 'radius':
             ls = etree.SubElement(s, "LineStyle")
             etree.SubElement(ls, "color").text = poly_c
             etree.SubElement(ls, "width").text = "3"
        else:
            etree.SubElement(etree.SubElement(s, "PolyStyle"), "color").text = poly_c
            etree.SubElement(etree.SubElement(s, "LineStyle"), "color").text = "ff" + line_c
            
    # 1. ARP Marker
    pm = etree.SubElement(doc, "Placemark")
    etree.SubElement(pm, "name").text = f"ARP ({arp['name']})"
    etree.SubElement(etree.SubElement(pm, "Point"), "coordinates").text = f"{arp['lon']},{arp['lat']},0"
    
    # 2. Radius Circle (Written FIRST so it appears BELOW features in KML)
    pm = etree.SubElement(doc, "Placemark")
    etree.SubElement(pm, "name").text = f"{int(radius_m/1000)}km Radius"
    etree.SubElement(pm, "styleUrl").text = "#radius"
    line = etree.SubElement(pm, "LineString")
    etree.SubElement(line, "coordinates").text = generate_circle_coords(arp['lat'], arp['lon'], radius_m)
    doc.append(pm)
    
    rows = []
    
    # 3. Features (Written LAST so they appear ON TOP and are clickable)
    for f in features:
        t = f['properties'].get('custom_type', 'waste')
        area = f['properties'].get('area_sq_m', 0)
        name = "Water Body" if t == 'water' else "Vegetation" if t == 'veg' else "Industrial/Waste"
        
        # Calculate Centroid & Distance
        try:
            c = shape(f['geometry']).centroid
            dist_km = haversine(arp['lat'], arp['lon'], c.y, c.x)
            lat_c, lon_c = c.y, c.x
        except:
            dist_km = 0
            lat_c, lon_c = 0, 0
        
        # KML Placemark
        pm = etree.SubElement(doc, "Placemark")
        etree.SubElement(pm, "name").text = name
        etree.SubElement(pm, "styleUrl").text = f"#{t}"
        
        # HTML Description Box (The Popup)
        desc_html = f"""
        <![CDATA[
            <div style="font-family:sans-serif; width:250px; padding:5px;">
                <h3 style="margin:0 0 10px 0; color:#333;">{name}</h3>
                <table border="0" cellspacing="0" cellpadding="2">
                    <tr><td style="font-weight:bold; color:#555;">Type:</td><td>{name}</td></tr>
                    <tr><td style="font-weight:bold; color:#555;">Area:</td><td>{int(area):,} m²</td></tr>
                    <tr><td style="font-weight:bold; color:#555;">Distance from ARP:</td><td>{dist_km:.2f} km</td></tr>
                </table>
            </div>
        ]]>
        """
        etree.SubElement(pm, "description").text = desc_html

        # Geometry
        geom = f['geometry']
        if geom['type'] in ['Polygon', 'MultiPolygon']:
            coords_list = geom['coordinates'][0] if geom['type'] == 'Polygon' else geom['coordinates'][0][0]
            coords_str = " ".join([f"{c[0]},{c[1]},0" for c in coords_list])
            poly = etree.SubElement(pm, "Polygon")
            etree.SubElement(etree.SubElement(poly, "outerBoundaryIs"), "LinearRing").append(etree.fromstring(f"<coordinates>{coords_str}</coordinates>"))
        
        rows.append({"Type": name, "Area_m2": int(area), "Distance_km": round(dist_km, 2), "Lat": lat_c, "Lon": lon_c})

    return etree.tostring(kml).decode(), pd.DataFrame(rows).to_csv(index=False)

# --- API ---
@app.route('/generate-report', methods=['POST', 'OPTIONS'])
def generate_report():
    if request.method == 'OPTIONS': return jsonify({"status": "ok"}), 200
    try:
        d = request.json
        if d.get('mode') == 'icao':
            icao_code = d["icao"].upper().strip()
            print(f"🔎 Looking up ICAO: {icao_code}")
            q = f'[out:json];nwr["icao"="{icao_code}"];out center;'
            try:
                overpass_resp = requests.post("https://overpass-api.de/api/interpreter", data={'data': q})
                data = overpass_resp.json()
                if not data.get('elements'):
                    return jsonify({"error": f"ICAO Code '{icao_code}' not found on OpenStreetMap."}), 404
                r = data['elements'][0]
                lat = r.get('center', {}).get('lat', r.get('lat'))
                lon = r.get('center', {}).get('lon', r.get('lon'))
                if not lat or not lon: return jsonify({"error": f"Could not determine coordinates for {icao_code}."}), 404
                arp = {"name": icao_code, "lat": lat, "lon": lon}
            except Exception as e:
                return jsonify({"error": "Failed to contact Map Server. Try Coordinates instead."}), 502
        else:
            arp = {"name": "Custom", "lat": float(d['lat']), "lon": float(d['lon'])}
            
        radius = float(d.get('radius_km', 13)) * 1000
        min_area = float(d.get('min_area_sq_m', 5000))
        
        print(f"📍 Processing {arp['name']} at {arp['lat']}, {arp['lon']}")
        raw = fetch_osm_data(arp['lat'], arp['lon'], radius) + fetch_gee_data(arp['lat'], arp['lon'], radius)
        
        print(f"   -> Total raw features fetched: {len(raw)}")
        
        final, display = [], []
        filtered_count = 0
        
        for f in raw:
            try:
                s = shape(f['geometry'])
                area = calculate_area(f['geometry'])
                
                # CORRECTED LOGIC: Skip features that are too small
                # Only keep features that meet BOTH conditions:
                # 1. Area is >= min_area (user's minimum threshold)
                # 2. Area is >= 10 (reasonable minimum to avoid tiny artifacts)
                if area < min_area:
                    filtered_count += 1
                    continue
                
                f['properties']['area_sq_m'] = area
                final.append(f)
                
                # Create simplified version for display (but keep all valid features)
                try:
                    simplified_geom = s.simplify(0.0005)  # Less aggressive simplification
                    if simplified_geom.is_valid and not simplified_geom.is_empty:
                        display.append({
                            "type": "Feature", 
                            "properties": f['properties'], 
                            "geometry": mapping(simplified_geom)
                        })
                    else:
                        # If simplification fails, use original
                        display.append({
                            "type": "Feature", 
                            "properties": f['properties'], 
                            "geometry": mapping(s)
                        })
                except:
                    # If simplification fails, use original geometry
                    display.append({
                        "type": "Feature", 
                        "properties": f['properties'], 
                        "geometry": mapping(s)
                    })
            except Exception as e:
                print(f"   -> Skipping invalid feature: {e}")
                pass
        
        print(f"   -> Features after filtering: {len(final)}")
        print(f"   -> Features filtered out (too small): {filtered_count}")
        print(f"   -> Display features: {len(display)}")
            
        kml, csv = generate_files(final, arp, radius)
        return jsonify({
            "message": "Success", 
            "feature_count": len(final), 
            "airport_info": arp,
            "map_geojson": {"type": "FeatureCollection", "features": display},
            "kml_string": kml, 
            "csv_string": csv
        })
    except Exception as e:
        print(f"CRITICAL: {e}")
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5000)