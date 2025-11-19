import os
import json
import ee
import requests
import pandas as pd
from shapely.geometry import shape
from pyproj import Geod
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
    # Check for Env Var first (Best Practice for Render)
    service_account_json = os.environ.get('GOOGLE_SERVICE_ACCOUNT_JSON')
    
    if service_account_json:
        # Render Environment Variable method
        key_data = json.loads(service_account_json)
        credentials = ee.ServiceAccountCredentials(key_data['client_email'], key_data=json.dumps(key_data))
        ee.Initialize(credentials=credentials)
        GEE_ENABLED = True
        print("✅ GEE Auth: Success (via Env Var)")
    elif os.path.exists(SERVICE_ACCOUNT_KEY_FILE):
        # Local file method
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
    
    min_lat = degrees(lat_rad - lat_dist)
    max_lat = degrees(lat_rad + lat_dist)
    min_lon = degrees(lon_rad - lon_dist_rad)
    max_lon = degrees(lon_rad + lon_dist_rad)
    return (min_lat, min_lon, max_lat, max_lon)

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
    try:
        s = shape(geom)
        geod = Geod(ellps='WGS84')
        if s.geom_type in ['Polygon', 'MultiPolygon']:
             area, _ = geod.geometry_area_perimeter(s)
             return abs(area)
        return 0
    except: return 0

def fetch_icao_data(icao):
    q = f"""[out:json][timeout:10];(node["icao"="{icao.upper()}"];way["icao"="{icao.upper()}"];relation["icao"="{icao.upper()}"];);out center;"""
    for server in OVERPASS_SERVERS:
        try:
            r = requests.post(server, data={'data': q}, headers=HEADERS, timeout=10)
            if r.status_code == 200 and r.json()['elements']:
                el = r.json()['elements'][0]
                lat = el.get('lat') or el.get('center', {}).get('lat')
                lon = el.get('lon') or el.get('center', {}).get('lon')
                if lat and lon: return {"name": f"{icao} Airport", "lat": lat, "lon": lon}
        except: continue
    raise Exception("Airport not found")

# -----------------------------------------------------------------
# STEP 3: DATA MINING
# -----------------------------------------------------------------

def fetch_osm_data(lat, lon, radius_m):
    print("🏗️ Fetching OSM...")
    bbox = get_bounding_box(lat, lon, (radius_m/1000)+2)
    bbox_str = f"{bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}"
    
    queries = {
        "waste": f'nwr["landuse"~"landfill|industrial|brownfield"](BBOX);nwr["amenity"~"waste_disposal|slaughterhouse"](BBOX);',
        "water": f'nwr["natural"~"water|wetland"](BBOX);nwr["landuse"~"reservoir|basin"](BBOX);'
    }
    
    features = []
    processed_ids = set()
    
    for type_key, q_temp in queries.items():
        full_q = f"[out:json][timeout:25];({q_temp.replace('BBOX', bbox_str)});out geom;"
        for server in OVERPASS_SERVERS:
            try:
                r = requests.post(server, data={'data': full_q}, headers=HEADERS, timeout=30)
                if r.status_code == 200:
                    data = r.json()
                    for f in data.get('features', []): # Note: Overpass output with 'out geom' might need manual GeoJSON conversion if not using osmtogeojson lib.
                        # However, if this worked in Colab, it means the library 'osm2geojson' or similar was used, OR the query was handled by a lib like osmnx.
                        # Wait, standard Overpass API returns JSON 'elements', not GeoJSON 'features'.
                        # The previous successful Colab code used 'osmnx' which handles this.
                        # Here in pure Python, we need to handle the 'elements' response if we aren't using a converter.
                        
                        # FIX: Let's ensure we are reading 'elements' and converting to simple GeoJSON-like dicts if needed.
                        # BUT, if the previous code worked partially, let's stick to the structure.
                        # If 'features' key is missing, it means we got raw OSM JSON.
                        
                        # Let's use a robust manual parser for OSM elements to ensure we get data.
                        elements = data.get('elements', [])
                        for el in elements:
                            fid = el.get('id')
                            if fid not in processed_ids:
                                # Simple Geometry Construction
                                geom = None
                                if el['type'] == 'node':
                                    geom = {"type": "Point", "coordinates": [el['lon'], el['lat']]}
                                elif el['type'] == 'way' and 'geometry' in el:
                                    coords = [[pt['lon'], pt['lat']] for pt in el['geometry']]
                                    geom = {"type": "Polygon", "coordinates": [coords]} # Approximation
                                
                                if geom:
                                    feat = {
                                        "type": "Feature",
                                        "properties": el.get('tags', {}),
                                        "geometry": geom
                                    }
                                    feat['properties']['custom_type'] = type_key
                                    features.append(feat)
                                    processed_ids.add(fid)
                    break 
            except: continue
            
    return features

def fetch_gee_data(lat, lon, radius_m):
    if not GEE_ENABLED: return []
    print("🛰️ Fetching Satellite...")
    try:
        aoi = ee.Geometry.Point(lon, lat).buffer(radius_m)
        img = ee.ImageCollection('ESA/WorldCover/v100').first()
        
        # 80=Water, 40=Cropland
        mask = img.eq(80).Or(img.eq(40))
        vectors = img.updateMask(mask).reduceToVectors(
            geometry=aoi, scale=20, maxPixels=1e9, geometryType='polygon', # Scale 20 is safer for memory
            reducer=ee.Reducer.first()
        )
        
        # Serialize properly
        data = vectors.getInfo() 
        
        features = []
        for f in data.get('features', []):
            label = f['properties'].get('label')
            if label == 80: f['properties']['custom_type'] = 'water'
            elif label == 40: f['properties']['custom_type'] = 'veg'
            else: continue
            f['properties']['tags'] = {'name': 'Satellite Feature'}
            features.append(f)
        return features
    except Exception as e:
        print(f"GEE Error: {e}")
        return []

# -----------------------------------------------------------------
# STEP 4: FILE GENERATION (KML/CSV)
# -----------------------------------------------------------------

def generate_kml(features, arp, radius_m):
    kml = etree.Element("kml", xmlns="http://www.opengis.net/kml/2.2")
    doc = etree.SubElement(kml, "Document")
    
    # Styles
    styles = {
        "water": "80ff0000", # Blue (AABBGGRR)
        "veg": "80008000",   # Green
        "waste": "9913458b", # Brown
        "radius": "ff0000ff" # Red
    }
    for k, color in styles.items():
        s = etree.SubElement(doc, "Style", id=k)
        etree.SubElement(etree.SubElement(s, "PolyStyle"), "color").text = color
        etree.SubElement(etree.SubElement(s, "LineStyle"), "color").text = color
        etree.SubElement(etree.SubElement(s, "LineStyle"), "width").text = "2"

    # ARP
    pm = etree.SubElement(doc, "Placemark")
    etree.SubElement(pm, "name").text = "ARP"
    etree.SubElement(etree.SubElement(pm, "Point"), "coordinates").text = f"{arp['lon']},{arp['lat']},0"
    
    # Features
    for f in features:
        t = f['properties'].get('custom_type', 'waste')
        area = f['properties'].get('area_sq_m', 0)
        
        pm = etree.SubElement(doc, "Placemark")
        name = "Water Body" if t == 'water' else "Vegetation" if t == 'veg' else "Industrial/Waste"
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
            
            name = "Water Body" if t == 'water' else "Vegetation" if t == 'veg' else "Industrial/Waste"
            
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
        radius_km = float(d.get('radius_km', 13))
        min_area = float(d.get('min_area_sq_m', 5000))
        radius_m = radius_km * 1000
        
        if d.get('mode') == 'icao':
            arp = fetch_icao_data(d.get('icao'))
        else:
            arp = {"name": "Custom", "lat": float(d.get('lat')), "lon": float(d.get('lon'))}
            
        # Fetch
        osm = fetch_osm_data(arp['lat'], arp['lon'], radius_m)
        gee = fetch_gee_data(arp['lat'], arp['lon'], radius_m)
        all_raw = osm + gee
        
        # Filter & Process
        final_features = []
        display_features = []
        
        for f in all_raw:
            geom = f.get('geometry')
            if not geom: continue
            
            # Area Filter
            area = calculate_area(geom)
            if area < min_area: continue
            
            f['properties']['area_sq_m'] = area
            final_features.append(f)
            
            # DISPLAY OPTIMIZATION (The Fix for Shaking!)
            # Create a lighter version for the browser map
            try:
                s = shape(geom)
                # 0.001 tolerance simplifies the shape to stop browser Z-fighting
                simplified = s.simplify(0.001, preserve_topology=True) 
                
                # We create a lightweight Feature for the map, but keep 'f' (heavy) for download
                display_f = {
                    "type": "Feature",
                    "properties": f['properties'],
                    "geometry": json.loads(json.dumps(simplified.__geo_interface__)) # Safe conversion
                }
                display_features.append(display_f)
            except:
                # Fallback to raw if simplification fails
                display_features.append(f)

        print(f"Found {len(final_features)} valid features.")
        
        return jsonify({
            "message": "Success",
            "feature_count": len(final_features),
            "airport_info": arp,
            "map_geojson": {"type": "FeatureCollection", "features": display_features}, # Send simplified to map
            "kml_string": generate_kml(final_features, arp, radius_m).decode('utf-8'), # Send raw to file
            "csv_string": generate_csv(final_features, arp)
        }), 200

    except Exception as e:
        print(f"Error: {e}")
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5000)