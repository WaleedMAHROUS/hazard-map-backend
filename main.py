import os
import json
import ee
import requests
import pandas as pd
from shapely.geometry import shape, mapping
from pyproj import Geod
from math import radians, cos, sin, asin, sqrt, degrees
from flask import Flask, request, jsonify
from flask_cors import CORS
from lxml import etree
import traceback

app = Flask(__name__)
CORS(app)

# --- AUTHENTICATION ---
GEE_ENABLED = False
try:
    # 1. Try Environment Variable (Best for Render)
    service_account_json = os.environ.get('GOOGLE_SERVICE_ACCOUNT_JSON')
    if service_account_json:
        print("🔑 Found GOOGLE_SERVICE_ACCOUNT_JSON env var.")
        key_data = json.loads(service_account_json)
        creds = ee.ServiceAccountCredentials(key_data['client_email'], key_data=json.dumps(key_data))
        ee.Initialize(credentials=creds)
        GEE_ENABLED = True
        print("✅ GEE Auth Success (Env Var)")
    
    # 2. Try Local File (Best for Local Testing)
    elif os.path.exists('service-account-key.json'):
        print("geefile found")
        creds = ee.ServiceAccountCredentials(json.load(open('service-account-key.json'))['client_email'], 'service-account-key.json')
        ee.Initialize(credentials=creds)
        GEE_ENABLED = True
        print("✅ GEE Auth Success (File)")
        
    else:
        print("⚠️ No Google Credentials found. Satellite data will be skipped.")
        
except Exception as e:
    print(f"❌ GEE Auth Failed: {e}")
    traceback.print_exc()

# --- HELPER FUNCTIONS ---
HEADERS = {'User-Agent': 'HabitatScanner/2.0'}
OVERPASS_SERVERS = [
    "https://overpass-api.de/api/interpreter",
    "https://lz4.overpass-api.de/api/interpreter",
    "https://z.overpass-api.de/api/interpreter"
]

def get_bbox(lat, lon, radius_km):
    R = 6371
    lat_rad, lon_rad = radians(lat), radians(lon)
    lat_d = radius_km / R
    lon_d = radius_km / (R * cos(lat_rad))
    return (degrees(lat_rad - lat_d), degrees(lon_rad - lon_d), degrees(lat_rad + lat_d), degrees(lon_rad + lon_d))

def haversine(lat1, lon1, lat2, lon2):
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

def fetch_icao_data(icao):
    print(f"🔎 Searching for ICAO: {icao}")
    q = f"""[out:json][timeout:15];(node["icao"="{icao.upper()}"];way["icao"="{icao.upper()}"];relation["icao"="{icao.upper()}"];);out center;"""
    
    for server in OVERPASS_SERVERS:
        try:
            print(f"   -> Trying server: {server}")
            r = requests.post(server, data={'data': q}, headers=HEADERS, timeout=15)
            if r.status_code == 200:
                data = r.json()
                if data.get('elements'):
                    el = data['elements'][0]
                    lat = el.get('lat') or el.get('center', {}).get('lat')
                    lon = el.get('lon') or el.get('center', {}).get('lon')
                    if lat and lon: 
                        print(f"✅ Found {icao} at {lat}, {lon}")
                        return {"name": f"{icao} Airport", "lat": lat, "lon": lon}
                else:
                     print(f"   -> No data found on {server}")
            else:
                print(f"   -> Server returned status {r.status_code}")
        except Exception as e: 
            print(f"   -> Error on {server}: {e}")
            continue
            
    # If loop finishes without returning, it failed everywhere
    raise Exception(f"Airport {icao} not found on any OSM server. Please check the code or use coordinates.")

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
    print("🏗️ Fetching OSM...")
    s, w, n, e = get_bbox(lat, lon, (radius_m/1000)+1)
    q = f"""[out:json][timeout:45];(nwr["landuse"~"landfill|industrial|brownfield|reservoir|basin"]({s},{w},{n},{e});nwr["amenity"~"waste_disposal|slaughterhouse"]({s},{w},{n},{e});nwr["natural"~"water|wetland"]({s},{w},{n},{e}););out body;>;out skel qt;"""
    
    for server in OVERPASS_SERVERS:
        try:
            resp = requests.post(server, data={'data': q}, headers=HEADERS, timeout=45)
            if resp.status_code != 200: continue
            
            feats = osm_to_geojson_custom(resp.json())
            for f in feats:
                p = f['properties']
                f['properties']['custom_type'] = 'water' if 'water' in str(p) or 'reservoir' in str(p) else 'waste'
                if f['geometry']['type'] == 'Point':
                    f['geometry'] = mapping(shape(f['geometry']).buffer(0.0005)) # Inflate points
            
            print(f"   -> OSM found {len(feats)}")
            return feats # Return immediately on success
        except Exception as e:
            print(f"OSM Query Error on {server}: {e}")
            continue
            
    print("❌ All OSM servers failed.")
    return []

def fetch_gee_data(lat, lon, radius_m):
    if not GEE_ENABLED: return []
    print("🛰️ Fetching Satellite...")
    try:
        aoi = ee.Geometry.Point(lon, lat).buffer(radius_m)
        img = ee.ImageCollection('ESA/WorldCover/v100').first()
        mask = img.eq(80).Or(img.eq(40))
        vectors = img.updateMask(mask).reduceToVectors(geometry=aoi, scale=20, maxPixels=1e9)
        feats = []
        for f in vectors.getInfo()['features']:
            f['properties']['custom_type'] = 'water' if f['properties']['label'] == 80 else 'veg'
            feats.append(f)
        print(f"   -> GEE found {len(feats)}")
        return feats
    except Exception as e:
        print(f"GEE Error: {e}"); return []

# --- FILES ---
def generate_files(features, arp):
    kml = etree.Element("kml", xmlns="http://www.opengis.net/kml/2.2")
    doc = etree.SubElement(kml, "Document")
    
    # Styles
    styles = {"water": "80ff0000", "veg": "80008000", "waste": "9913458b"}
    for k, c in styles.items():
        s = etree.SubElement(doc, "Style", id=k)
        etree.SubElement(etree.SubElement(s, "PolyStyle"), "color").text = c
    
    rows = []
    for f in features:
        t = f['properties'].get('custom_type', 'waste')
        area = f['properties'].get('area_sq_m', 0)
        name = "Water" if t == 'water' else "Veg" if t == 'veg' else "Waste"
        
        pm = etree.SubElement(doc, "Placemark")
        etree.SubElement(pm, "name").text = name
        etree.SubElement(pm, "styleUrl").text = f"#{t}"
        
        # CSV
        try:
            c = shape(f['geometry']).centroid
            rows.append({"Type": name, "Area": area, "Lat": c.y, "Lon": c.x})
        except: pass
        
    return etree.tostring(kml).decode(), pd.DataFrame(rows).to_csv(index=False)

# --- API ---
@app.route('/generate-report', methods=['POST', 'OPTIONS'])
def generate_report():
    if request.method == 'OPTIONS': return jsonify({"status": "ok"}), 200
    try:
        d = request.json
        radius = float(d.get('radius_km', 13)) * 1000
        min_area = float(d.get('min_area_sq_m', 5000))

        if d.get('mode') == 'icao':
            arp = fetch_icao_data(d.get('icao'))
        else:
            arp = {"name": "Custom", "lat": float(d['lat']), "lon": float(d['lon'])}
            
        raw = fetch_osm_data(arp['lat'], arp['lon'], radius) + fetch_gee_data(arp['lat'], arp['lon'], radius)
        
        final, display = [], []
        for f in raw:
            s = shape(f['geometry'])
            area = calculate_area(f['geometry'])
            if area < min_area and area > 10: continue 
            f['properties']['area_sq_m'] = area
            final.append(f)
            display.append({"type": "Feature", "properties": f['properties'], "geometry": mapping(s.simplify(0.001))})
            
        kml, csv = generate_files(final, arp)
        return jsonify({
            "message": "Success", "feature_count": len(final), "airport_info": arp,
            "map_geojson": {"type": "FeatureCollection", "features": display},
            "kml_string": kml, "csv_string": csv
        })
    except Exception as e:
        print(f"CRITICAL: {e}")
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5000)