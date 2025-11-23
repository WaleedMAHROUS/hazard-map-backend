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
from flask_caching import Cache
from lxml import etree
import traceback

app = Flask(__name__)
CORS(app)
cache = Cache(app, config={'CACHE_TYPE': 'SimpleCache', 'CACHE_DEFAULT_TIMEOUT': 3600})

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

def calculate_risk_score(feature, arp_lat, arp_lon):
    """
    Calculates a risk score (1-10) based on:
    1. Distance from ARP (Closer = Higher Risk)
    2. Area Size (Larger = Higher Risk)
    3. Feature Type (Specific types = Higher Risk)
    """
    try:
        # 1. Distance Score (0-4 points)
        # < 3km: 4 pts, < 6km: 3 pts, < 9km: 2 pts, < 13km: 1 pt
        c = shape(feature['geometry']).centroid
        dist_km = haversine(arp_lat, arp_lon, c.y, c.x)
        
        dist_score = 0
        if dist_km < 3: dist_score = 4
        elif dist_km < 6: dist_score = 3
        elif dist_km < 9: dist_score = 2
        elif dist_km < 13: dist_score = 1
        
        # 2. Area Score (0-3 points)
        # > 1km²: 3 pts, > 0.1km²: 2 pts, > 0.01km²: 1 pt
        area_sq_m = feature['properties'].get('area_sq_m', 0)
        area_score = 0
        if area_sq_m > 1000000: area_score = 3
        elif area_sq_m > 100000: area_score = 2
        elif area_sq_m > 10000: area_score = 1
        
        # 3. Type Score (0-3 points)
        # Waste/Landfill: 3 pts, Water: 2 pts, Veg: 1 pt
        type_score = 1 # Default (Veg/Other)
        t = feature['properties'].get('custom_type', '')
        if t == 'waste': type_score = 3
        elif t == 'water': type_score = 2
        
        total_score = dist_score + area_score + type_score
        return min(10, max(1, total_score)) # Clamp between 1 and 10
    except:
        return 1

def calculate_distance(feature, arp_lat, arp_lon):
    try:
        c = shape(feature['geometry']).centroid
        return haversine(arp_lat, arp_lon, c.y, c.x)
    except:
        return 0

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
    # Refined Query for HIGH-RISK wildlife attractants only
    q = f"""[out:json][timeout:90];
    (
      nwr["landuse"~"landfill|industrial|brownfield|reservoir|basin|farmland|orchard|vineyard"]["aeroway"!~"."]({s},{w},{n},{e});
      nwr["amenity"~"waste_disposal|slaughterhouse"]["aeroway"!~"."]({s},{w},{n},{e});
      nwr["natural"~"water|wetland|scrub|heath"]["aeroway"!~"."]({s},{w},{n},{e});
    );
    out body;
    >;
    out skel qt;"""
    try:
        resp = requests.post("https://overpass-api.de/api/interpreter", data={'data': q}, headers=HEADERS)
        if resp.status_code != 200: return []
        feats = osm_to_geojson_custom(resp.json())
        for f in feats:
            p = f['properties']
            
            # Proper classification based on OSM tags
            landuse = p.get('landuse', '').lower()
            natural = p.get('natural', '').lower()
            leisure = p.get('leisure', '').lower()
            amenity = p.get('amenity', '').lower()
            aeroway = p.get('aeroway', '').lower()
            
            # Skip airport infrastructure
            if aeroway in ['runway', 'taxiway', 'apron', 'terminal', 'hangar', 'gate']:
                continue
            
            # Classify as waste/industrial
            if landuse in ['landfill', 'industrial', 'brownfield'] or amenity in ['waste_disposal', 'slaughterhouse']:
                f['properties']['custom_type'] = 'waste'
            # Classify as water
            elif 'water' in natural or 'water' in landuse or 'reservoir' in landuse or 'wetland' in natural or 'basin' in landuse:
                f['properties']['custom_type'] = 'water'
            # Classify as vegetation (everything else: parks, farms, grassland, etc.)
            else:
                f['properties']['custom_type'] = 'veg'
            
            f['properties']['source'] = 'OpenStreetMap'
            if f['geometry']['type'] == 'Point':
                # FIX: Reduced buffer from 0.0005 (~9600m2) to 0.0001 (~380m2)
                # This ensures points aren't artificially larger than the 5000m2 filter
                f['geometry'] = mapping(shape(f['geometry']).buffer(0.0001)) 
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
        
        # OPTIMIZED: Only detect water (80), not vegetation (40)
        # This makes satellite data complement OSM rather than overwhelm it
        mask = img.eq(80)
        
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
            # All GEE features are now water (since we only masked for class 80)
            f['properties']['custom_type'] = 'water'
            f['properties']['source'] = 'Satellite (GEE)'
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
            etree.SubElement(etree.SubElement(s, "LineStyle"), "width").text = "2" # Thicker lines for visibility
            
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
        risk = f['properties'].get('risk_score', 1)
        dist_km = f['properties'].get('dist_km', 0)
        
        # Use specific name if available, else generic
        specific_name = f['properties'].get('name', '')
        generic_name = "Water Body" if t == 'water' else "Vegetation" if t == 'veg' else "Industrial/Waste"
        display_name = specific_name if specific_name else generic_name
        
        # KML Placemark
        pm = etree.SubElement(doc, "Placemark")
        etree.SubElement(pm, "name").text = f"{display_name} (Risk: {risk})"
        etree.SubElement(pm, "styleUrl").text = f"#{t}"
        
        # HTML Description Box (The Popup)
        desc_html = f"""
        <![CDATA[
            <div style="font-family:sans-serif; width:300px; padding:10px;">
                <h3 style="margin:0 0 10px 0; color:#333; border-bottom:1px solid #ccc; padding-bottom:5px;">{display_name}</h3>
                <table border="0" cellspacing="0" cellpadding="4" style="font-size:12px;">
                    <tr><td style="font-weight:bold; color:#555;">Type:</td><td>{generic_name}</td></tr>
                    <tr><td style="font-weight:bold; color:#555;">Source:</td><td>{f['properties'].get('source', 'Unknown')}</td></tr>
                    <tr><td style="font-weight:bold; color:#555;">Risk Score:</td><td><b style="color:{'red' if risk >= 7 else 'orange' if risk >= 4 else 'green'};">{risk}/10</b></td></tr>
                    <tr><td style="font-weight:bold; color:#555;">Area:</td><td>{int(area):,} m²</td></tr>
                    <tr><td style="font-weight:bold; color:#555;">Distance from ARP:</td><td>{dist_km:.2f} km</td></tr>
                </table>
                <p style="font-size:10px; color:#777; margin-top:10px;">Detected by Habitat Scanner v0.6</p>
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
        
        # CSV Row - Renamed Hazards to Habitats
        rows.append({
            "Habitat Name": display_name,
            "Habitat Type": generic_name, 
            "Source": f['properties'].get('source', 'Unknown'),
            "Risk Score": risk,
            "Area (m2)": int(area), 
            "Distance (km)": round(dist_km, 2), 
            "Lat": f['properties'].get('lat_c', 0), 
            "Lon": f['properties'].get('lon_c', 0)
        })

    return etree.tostring(kml).decode(), pd.DataFrame(rows).to_csv(index=False)

# --- API ---
@app.route('/generate-report', methods=['POST', 'OPTIONS'])
# @cache.cached(timeout=3600, key_prefix=lambda: request.data) # Cache based on request body
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
        min_area_input = d.get('min_area_sq_m', 0)
        
        # Handle empty string or None
        if min_area_input == '' or min_area_input is None:
            min_area = 0
        else:
            min_area = float(min_area_input)
        
        print(f"📍 Processing {arp['name']} at {arp['lat']}, {arp['lon']}")
        print(f"   -> Minimum area threshold: {min_area} m² (type: {type(min_area).__name__})")
        osm_data = fetch_osm_data(arp['lat'], arp['lon'], radius)
        gee_data = fetch_gee_data(arp['lat'], arp['lon'], radius)
        raw = osm_data + gee_data
        
        print(f"   -> Total raw features fetched: {len(raw)}")
        
        final, display = [], []
        filtered_count = 0
        invalid_count = 0
        
        for f in raw:
            try:
                s = shape(f['geometry'])
                
                # Attempt to repair invalid geometry using buffer(0) trick
                if not s.is_valid:
                    try:
                        s = s.buffer(0)  # This fixes most topology issues
                        if s.is_valid and not s.is_empty:
                            f['geometry'] = mapping(s)
                        else:
                            invalid_count += 1
                            continue
                    except:
                        invalid_count += 1
                        continue
                
                # Skip empty geometries
                if s.is_empty:
                    invalid_count += 1
                    continue
                
                area = calculate_area(f['geometry'])
                
                # Debug: Log first few features
                if len(final) < 3:
                    print(f"   -> Feature area: {area:.2f} m², min_area: {min_area}, will_filter: {min_area > 0 and area < min_area}")
                
                # FIXED LOGIC: Only filter based on min_area threshold
                # If min_area is 0, keep all valid features regardless of size
                if min_area > 0 and area < min_area:
                    filtered_count += 1
                    continue
                
                # Store the calculated area
                f['properties']['area_sq_m'] = area
                
                # Calculate Risk Score
                risk_score = calculate_risk_score(f, arp['lat'], arp['lon'])
                f['properties']['risk_score'] = risk_score

                # Calculate Distance
                dist_km = calculate_distance(f, arp['lat'], arp['lon'])
                f['properties']['dist_km'] = dist_km
                
                # Store Centroid for CSV
                try:
                    c = shape(f['geometry']).centroid
                    f['properties']['lat_c'] = c.y
                    f['properties']['lon_c'] = c.x
                except: pass
                
                final.append(f)
                
                # Create simplified version for display
                try:
                    # Use less aggressive simplification to preserve features
                    simplified_geom = s.simplify(0.0003, preserve_topology=True)
                    if simplified_geom.is_valid and not simplified_geom.is_empty:
                        display.append({
                            "type": "Feature", 
                            "properties": f['properties'], 
                            "geometry": mapping(simplified_geom)
                        })
                    else:
                        # If simplification creates invalid geometry, use original
                        display.append({
                            "type": "Feature", 
                            "properties": f['properties'], 
                            "geometry": mapping(s)
                        })
                except Exception as simp_err:
                    # If simplification fails, use original geometry
                    print(f"   -> Simplification failed, using original: {simp_err}")
                    display.append({
                        "type": "Feature", 
                        "properties": f['properties'], 
                        "geometry": mapping(s)
                    })
            except Exception as e:
                invalid_count += 1
                print(f"   -> Skipping invalid feature: {e}")
                pass
        
        print(f"   -> Features after filtering: {len(final)}")
        print(f"   -> Features filtered (too small): {filtered_count}")
        print(f"   -> Features invalid (bad geometry): {invalid_count}")
        print(f"   -> Display features: {len(display)}")
            
        kml, csv = generate_files(final, arp, radius)
        
        # Calculate Statistics
        stats = {
            "total_area": sum(f['properties']['area_sq_m'] for f in final),
            "by_type": {},
            "by_risk": {"High": 0, "Medium": 0, "Low": 0}
        }
        
        for f in final:
            t = f['properties'].get('custom_type', 'other')
            stats["by_type"][t] = stats["by_type"].get(t, 0) + 1
            
            rs = f['properties'].get('risk_score', 1)
            if rs >= 7: stats["by_risk"]["High"] += 1
            elif rs >= 4: stats["by_risk"]["Medium"] += 1
            else: stats["by_risk"]["Low"] += 1

        return jsonify({
            "message": "Success", 
            "feature_count": len(final), 
            "counts": {
                "osm": len(osm_data),
                "gee": len(gee_data),
                "filtered": filtered_count,
                "invalid": invalid_count,
                "final": len(final)
            },
            "airport_info": arp,
            "map_geojson": {"type": "FeatureCollection", "features": display},
            "stats": stats,
            "kml_string": kml, 
            "csv_string": csv
        })
    except Exception as e:
        print(f"CRITICAL: {e}")
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5000)