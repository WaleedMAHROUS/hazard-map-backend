import os
import json
import ee
import requests
import pandas as pd
from shapely.geometry import shape
from shapely.ops import transform
from pyproj import Geod, Transformer
from math import radians, cos, sin, asin, sqrt, degrees, atan2
from flask import Flask, request, jsonify
from flask_cors import CORS
from lxml import etree
import datetime

# -----------------------------------------------------------------
# STEP 1: SETUP & AUTHENTICATION
# -----------------------------------------------------------------

app = Flask(__name__)
CORS(app)

SERVICE_ACCOUNT_KEY_FILE = 'service-account-key.json'
GEE_ENABLED = False

try:
    service_account_json = os.environ.get('GOOGLE_SERVICE_ACCOUNT_JSON')
    if service_account_json:
        key_data = json.loads(service_account_json)
        credentials = ee.ServiceAccountCredentials(key_data['client_email'], key_data=json.dumps(key_data))
        ee.Initialize(credentials=credentials)
        GEE_ENABLED = True
        print("✅ GEE Auth: Success (via Env Var)")
    elif os.path.exists(SERVICE_ACCOUNT_KEY_FILE):
        with open(SERVICE_ACCOUNT_KEY_FILE, 'r') as f:
            key_data = json.load(f)
        credentials = ee.ServiceAccountCredentials(key_data['client_email'], SERVICE_ACCOUNT_KEY_FILE)
        ee.Initialize(credentials=credentials)
        GEE_ENABLED = True
        print("✅ GEE Auth: Success (via File)")
    else:
        print("⚠️ GEE Auth: No key found. Satellite data disabled.")
except Exception as e:
    print(f"❌ GEE Auth Failed: {e}")

# -----------------------------------------------------------------
# STEP 2: HELPER FUNCTIONS
# -----------------------------------------------------------------

OVERPASS_SERVERS = [
    "https://overpass-api.de/api/interpreter",
    "https://lz4.overpass-api.de/api/interpreter"
]
HEADERS = {'User-Agent': 'HabitatScanner/Pro'}

def get_bounding_box(lat, lon, radius_km):
    R = 6371
    lat_rad, lon_rad = radians(lat), radians(lon)
    lat_dist = radius_km / R
    lon_dist_rad = radius_km / (R * cos(lat_rad))
    return (
        degrees(lat_rad - lat_dist),
        degrees(lon_rad - lon_dist_rad),
        degrees(lat_rad + lat_dist),
        degrees(lon_rad + lon_dist_rad)
    )

def haversine_distance(lat1, lon1, lat2, lon2):
    R = 6371
    dLat, dLon = radians(lat2 - lat1), radians(lon2 - lon1)
    lat1, lat2 = radians(lat1), radians(lat2)
    a = sin(dLat / 2)**2 + cos(lat1) * cos(lat2) * sin(dLon / 2)**2
    c = 2 * asin(sqrt(a))
    return R * c

def get_centroid(geom):
    s = shape(geom)
    return s.centroid.y, s.centroid.x

def calculate_area(geom):
    """
    Robust area calculation using pyproj Geod (ellipsoidal).
    Returns area in square meters.
    """
    try:
        s = shape(geom)
        if s.is_empty: return 0
        
        # Geodetic calculation is most accurate for lat/lon data
        geod = Geod(ellps='WGS84')
        if s.geom_type in ['Polygon', 'MultiPolygon']:
             area, _ = geod.geometry_area_perimeter(s)
             return abs(area)
        return 0
    except Exception as e:
        print(f"⚠️ Area Calc Error: {e}")
        return 0

def fetch_icao_data(icao):
    print(f"🔎 Searching for ICAO: {icao}")
    q = f"""[out:json][timeout:10];(node["icao"="{icao.upper()}"];way["icao"="{icao.upper()}"];relation["icao"="{icao.upper()}"];);out center;"""
    for server in OVERPASS_SERVERS:
        try:
            r = requests.post(server, data={'data': q}, headers=HEADERS, timeout=10)
            if r.status_code == 200 and r.json()['elements']:
                el = r.json()['elements'][0]
                lat = el.get('lat') or el.get('center', {}).get('lat')
                lon = el.get('lon') or el.get('center', {}).get('lon')
                if lat and lon: 
                    print(f"✅ Found {icao} at {lat}, {lon}")
                    return {"name": f"{icao} Airport", "lat": lat, "lon": lon}
        except Exception as e: 
            print(f"Server Error: {e}")
            continue
    raise Exception("Airport not found")

# -----------------------------------------------------------------
# STEP 3: DATA MINING
# -----------------------------------------------------------------

def fetch_osm_data(lat, lon, radius_m):
    print("🏗️ Fetching OSM Data...")
    # Bounding box: (South, West, North, East) for Overpass
    # Note: get_bounding_box returns (min_lat, min_lon, max_lat, max_lon)
    min_lat, min_lon, max_lat, max_lon = get_bounding_box(lat, lon, (radius_m/1000)+2)
    bbox_str = f"{min_lat},{min_lon},{max_lat},{max_lon}"
    
    queries = {
        "waste": f'nwr["landuse"~"landfill|industrial|brownfield"]({bbox_str});nwr["amenity"~"waste_disposal|slaughterhouse"]({bbox_str});',
        "water": f'nwr["natural"~"water|wetland"]({bbox_str});nwr["landuse"~"reservoir|basin"]({bbox_str});'
    }
    
    features = []
    processed_ids = set()
    
    for type_key, q_temp in queries.items():
        full_q = f"[out:json][timeout:25];({q_temp});out geom;"
        for server in OVERPASS_SERVERS:
            try:
                r = requests.post(server, data={'data': full_q}, headers=HEADERS, timeout=30)
                if r.status_code == 200:
                    data = r.json()
                    count = 0
                    for f in data.get('features', []): # Note: Overpass usually returns 'elements', need conversion to GeoJSON?
                        # The python `requests` response from Overpass is NOT GeoJSON by default.
                        # It's a custom JSON format unless we use a library or specific output.
                        # Wait! Your previous successful code relied on Overpass returning GeoJSON?
                        # Overpass API does NOT return GeoJSON directly unless using `osmtogeojson` or similar client-side.
                        # BUT wait, `overpass-api.de/api/interpreter` returns { "elements": [...] }
                        # We need to convert this to GeoJSON features manually here.
                        pass # Logic handled below
                    
                    # Actually, let's use a cleaner method. 
                    # If we are raw querying, we get Elements.
                    # Let's trust the 'elements' structure and convert manually if needed.
                    # HOWEVER, to match your successful Colab run (which used osmnx), we should stick to logic that works.
                    # Since we can't easily install osmnx on standard Render without GDAL hell, 
                    # let's ensure we parse the Overpass JSON correctly.
                    
                    # FIX: We will use `osm2geojson` library if available, or manual parsing.
                    # Since we don't want to add complex deps, let's manually parse basic OSM elements to GeoJSON-like dicts.
                    
                    elements = data.get('elements', [])
                    print(f"   -> Raw OSM elements found for {type_key}: {len(elements)}")
                    
                    for el in elements:
                        fid = el.get('id')
                        if fid not in processed_ids:
                            # Simple geometry construction
                            geom = None
                            if el['type'] == 'node':
                                geom = {"type": "Point", "coordinates": [el['lon'], el['lat']]}
                            elif el['type'] == 'way' and 'geometry' in el:
                                # Ways with 'out geom' have a 'geometry' list of dicts {lat, lon}
                                coords = [[pt['lon'], pt['lat']] for pt in el['geometry']]
                                # Close polygon if it's a polygon type
                                if coords[0] == coords[-1]:
                                    geom = {"type": "Polygon", "coordinates": [coords]}
                                else:
                                    geom = {"type": "LineString", "coordinates": coords}
                            elif 'bounds' in el: # Relations sometimes
                                # Centroid fallback for complex relations
                                b = el['bounds']
                                c_lat = (b['minlat'] + b['maxlat']) / 2
                                c_lon = (b['minlon'] + b['maxlon']) / 2
                                geom = {"type": "Point", "coordinates": [c_lon, c_lat]}
                                
                            if geom:
                                feat = {
                                    "type": "Feature",
                                    "id": fid,
                                    "properties": el.get('tags', {}),
                                    "geometry": geom
                                }
                                feat['properties']['custom_type'] = type_key
                                features.append(feat)
                                processed_ids.add(fid)
                                count += 1
                    print(f"   -> Converted {count} {type_key} features.")
                    break 
            except Exception as e: 
                print(f"   -> OSM Query Error: {e}")
                continue
            
    return features

def fetch_gee_data(lat, lon, radius_m):
    if not GEE_ENABLED: return []
    print("🛰️ Fetching Satellite Data...")
    try:
        aoi = ee.Geometry.Point(lon, lat).buffer(radius_m)
        img = ee.ImageCollection('ESA/WorldCover/v100').first()
        
        mask = img.eq(80).Or(img.eq(40))
        vectors = img.updateMask(mask).reduceToVectors(
            geometry=aoi, scale=20, maxPixels=1e9, geometryType='polygon', # Scale 20 is safer for memory
            reducer=ee.Reducer.first()
        )
        
        # Need to serialize EE object to Python Dict
        data = vectors.getInfo() 
        
        features = []
        for f in data.get('features', []):
            label = f['properties'].get('label')
            if label == 80: f['properties']['custom_type'] = 'water'
            elif label == 40: f['properties']['custom_type'] = 'veg'
            else: continue
            f['properties']['tags'] = {'name': 'Satellite Feature'}
            features.append(f)
            
        print(f"   -> Found {len(features)} Satellite features.")
        return features
    except Exception as e:
        print(f"❌ GEE Error: {e}")
        return []

# -----------------------------------------------------------------
# STEP 4: FILE GENERATION
# -----------------------------------------------------------------

def generate_kml(features, arp, radius_m):
    kml = etree.Element("kml", xmlns="http://www.opengis.net/kml/2.2")
    doc = etree.SubElement(kml, "Document")
    
    styles = {"water": "80ff0000", "veg": "80008000", "waste": "9913458b", "radius": "ff0000ff"}
    for k, color in styles.items():
        s = etree.SubElement(doc, "Style", id=k)
        etree.SubElement(etree.SubElement(s, "PolyStyle"), "color").text = color
        etree.SubElement(etree.SubElement(s, "LineStyle"), "color").text = color
        etree.SubElement(etree.SubElement(s, "LineStyle"), "width").text = "2"

    pm = etree.SubElement(doc, "Placemark")
    etree.SubElement(pm, "name").text = "ARP"
    etree.SubElement(etree.SubElement(pm, "Point"), "coordinates").text = f"{arp['lon']},{arp['lat']},0"
    
    for f in features:
        t = f['properties'].get('custom_type', 'waste')
        area = f['properties'].get('area_sq_m', 0)
        pm = etree.SubElement(doc, "Placemark")
        name = "Water Body" if t == 'water' else "Vegetation" if t == 'veg' else "Ind./Waste"
        etree.SubElement(pm, "name").text = name
        etree.SubElement(pm, "description").text = f"Area: {int(area):,} m2"
        etree.SubElement(pm, "styleUrl").text = f"#{t}"
        
        geom = f['geometry']
        if geom['type'] == 'Polygon':
            coords = " ".join([f"{c[0]},{c[1]},0" for c in geom['coordinates'][0]])
            poly = etree.SubElement(pm, "Polygon")
            ob = etree.SubElement(poly, "outerBoundaryIs")
            lr = etree.SubElement(ob, "LinearRing")
            etree.SubElement(lr, "coordinates").text = coords

    return etree.tostring(kml, pretty_print=True)

def generate_csv(features, arp):
    rows = []
    for f in features:
        try:
            c_lat, c_lon = get_centroid(f['geometry'])
            dist = haversine_distance(arp['lat'], arp['lon'], c_lat, c_lon)
            t = f['properties'].get('custom_type')
            name = "Water Body" if t == 'water' else "Vegetation" if t == 'veg' else "Ind./Waste"
            rows.append({
                "Type": name,
                "Area_m2": int(f['properties']['area_sq_m']),
                "Distance_km": round(dist, 2),
                "Lat": c_lat, "Lon": c_lon
            })
        except: continue
    
    df = pd.DataFrame(rows)
    return df.to_csv(index=False) if not df.empty else ""

# -----------------------------------------------------------------
# STEP 5: ENDPOINT
# -----------------------------------------------------------------

@app.route('/generate-report', methods=['POST', 'OPTIONS'])
def generate_report():
    if request.method == 'OPTIONS': return jsonify({"status": "ok"}), 200
    
    try:
        d = request.json
        print(f"📥 Received Request: {d}")
        
        radius_km = float(d.get('radius_km', 13))
        min_area = float(d.get('min_area_sq_m', 5000))
        radius_m = radius_km * 1000
        
        if d.get('mode') == 'icao':
            arp = fetch_icao_data(d.get('icao'))
        else:
            arp = {"name": "Custom", "lat": float(d.get('lat')), "lon": float(d.get('lon'))}
            
        print(f"📍 Center: {arp}")

        # Fetch
        osm = fetch_osm_data(arp['lat'], arp['lon'], radius_m)
        gee = fetch_gee_data(arp['lat'], arp['lon'], radius_m)
        all_raw = osm + gee
        
        print(f"📊 Raw Features Found: {len(all_raw)}")
        
        # Filter & Process
        final_features = []
        display_features = []
        
        for f in all_raw:
            geom = f.get('geometry')
            if not geom: continue
            
            # Area Filter
            area = calculate_area(geom)
            # print(f"   - Feature Area: {area}") # Uncomment for deep debugging
            
            if area < min_area: continue
            
            f['properties']['area_sq_m'] = area
            final_features.append(f)
            
            # DISPLAY OPTIMIZATION
            try:
                s = shape(geom)
                simplified = s.simplify(0.001, preserve_topology=True) 
                display_f = {
                    "type": "Feature",
                    "properties": f['properties'],
                    "geometry": json.loads(json.dumps(simplified.__geo_interface__))
                }
                display_features.append(display_f)
            except:
                display_features.append(f)

        print(f"✅ Final Filtered Features: {len(final_features)}")
        
        return jsonify({
            "message": "Success",
            "feature_count": len(final_features),
            "airport_info": arp,
            "map_geojson": {"type": "FeatureCollection", "features": display_features},
            "kml_string": generate_kml(final_features, arp, radius_m).decode('utf-8'),
            "csv_string": generate_csv(final_features, arp)
        }), 200

    except Exception as e:
        print(f"🔥 CRITICAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5000)