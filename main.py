import os
import json
import ee
import requests
import pandas as pd
from shapely.geometry import shape, mapping, Point, LineString, Polygon
from pyproj import Geod
from math import radians, cos, sin, asin, sqrt, degrees
from flask import Flask, request, jsonify
from flask_cors import CORS
from lxml import etree

app = Flask(__name__)
CORS(app)

# --- AUTHENTICATION ---
GEE_ENABLED = False
try:
    service_account_json = os.environ.get('GOOGLE_SERVICE_ACCOUNT_JSON')
    if service_account_json:
        creds = ee.ServiceAccountCredentials(json.loads(service_account_json)['client_email'], key_data=service_account_json)
        ee.Initialize(credentials=creds)
        GEE_ENABLED = True
        print("✅ GEE Auth: Success (Env Var)")
    elif os.path.exists('service-account-key.json'):
        creds = ee.ServiceAccountCredentials(json.load(open('service-account-key.json'))['client_email'], 'service-account-key.json')
        ee.Initialize(credentials=creds)
        GEE_ENABLED = True
        print("✅ GEE Auth: Success (File)")
except Exception as e:
    print(f"⚠️ GEE Auth Failed: {e}")

# --- HELPER FUNCTIONS ---
HEADERS = {'User-Agent': 'HabitatScanner/2.0'}

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

# --- ROBUST CUSTOM PARSER (Replaces osm2geojson) ---
def osm_to_geojson_custom(osm_json):
    """
    Converts raw OSM JSON (elements) to GeoJSON Features without external libs.
    Handles Nodes, Ways, and Relations (Multipolygons).
    """
    elements = osm_json.get('elements', [])
    nodes = {e['id']: (e['lon'], e['lat']) for e in elements if e['type'] == 'node'}
    ways = {e['id']: e for e in elements if e['type'] == 'way'}
    
    features = []
    
    for el in elements:
        # 1. STANDALONE NODES (Points)
        if el['type'] == 'node' and 'tags' in el:
            feat = {
                "type": "Feature",
                "properties": el['tags'],
                "geometry": {"type": "Point", "coordinates": [el['lon'], el['lat']]}
            }
            features.append(feat)
            
        # 2. WAYS (Polygons or Lines)
        elif el['type'] == 'way' and 'tags' in el:
            coords = []
            for nid in el.get('nodes', []):
                if nid in nodes: coords.append(nodes[nid])
                
            if len(coords) > 2:
                # Check if closed (Polygon)
                if coords[0] == coords[-1]:
                    geom_type = "Polygon"
                    geom_coords = [coords]
                else:
                    geom_type = "LineString"
                    geom_coords = coords
                    
                feat = {
                    "type": "Feature",
                    "properties": el['tags'],
                    "geometry": {"type": geom_type, "coordinates": geom_coords}
                }
                features.append(feat)
                
        # 3. RELATIONS (Multipolygons) - Simplified
        elif el['type'] == 'relation' and el.get('tags', {}).get('type') == 'multipolygon':
            # Basic relation handler: Grab "outer" ways and merge them
            # Complex multipolygon stitching is hard in pure python, 
            # so we fallback to taking the first valid "outer" way found as the representation.
            # This covers 90% of standard use cases like airports/parks.
            outer_coords = []
            for member in el.get('members', []):
                if member['role'] == 'outer' and member['type'] == 'way':
                    wid = member['ref']
                    if wid in ways:
                        w = ways[wid]
                        w_coords = []
                        for nid in w.get('nodes', []):
                            if nid in nodes: w_coords.append(nodes[nid])
                        if w_coords:
                            outer_coords = w_coords
                            break # Take first valid outer ring
            
            if outer_coords:
                 feat = {
                    "type": "Feature",
                    "properties": el['tags'],
                    "geometry": {"type": "Polygon", "coordinates": [outer_coords]}
                }
                 features.append(feat)

    return features

# --- DATA MINING ---
def fetch_osm_data(lat, lon, radius_m):
    print("🏗️ Fetching OSM...")
    s, w, n, e = get_bbox(lat, lon, (radius_m/1000)+1)
    bbox = f"{s},{w},{n},{e}"
    
    # Recursive query to get nodes/ways for geometry construction
    q = f"""
    [out:json][timeout:45];
    (
      nwr["landuse"~"landfill|industrial|brownfield|reservoir|basin"]({bbox});
      nwr["amenity"~"waste_disposal|slaughterhouse"]({bbox});
      nwr["natural"~"water|wetland"]({bbox});
    );
    out body;
    >;
    out skel qt;
    """
    
    try:
        resp = requests.post("https://overpass-api.de/api/interpreter", data={'data': q}, headers=HEADERS)
        if resp.status_code != 200: return []
        
        # Use our custom robust parser
        raw_features = osm_to_geojson_custom(resp.json())
        
        features = []
        for f in raw_features:
            props = f.get('properties', {})
            
            # Determine Type
            f_type = 'waste'
            if 'water' in props.get('natural', '') or 'water' in props.get('landuse', ''): f_type = 'water'
            elif 'reservoir' in props.get('landuse', ''): f_type = 'water'
            elif 'wetland' in props.get('natural', ''): f_type = 'water'
            
            f['properties']['custom_type'] = f_type
            
            # Point Inflation: Buffer Points to 50m Circles so they pass area filter
            if f['geometry']['type'] == 'Point':
                s = shape(f['geometry'])
                # Buffer approx 50m (0.00045 deg)
                f['geometry'] = mapping(s.buffer(0.00045))
                
            features.append(f)
            
        print(f"   -> OSM found {len(features)} items")
        return features
    except Exception as e:
        print(f"OSM Error: {e}")
        return []

def fetch_gee_data(lat, lon, radius_m):
    if not GEE_ENABLED: return []
    print("🛰️ Fetching Satellite...")
    try:
        aoi = ee.Geometry.Point(lon, lat).buffer(radius_m)
        img = ee.ImageCollection('ESA/WorldCover/v100').first()
        mask = img.eq(80).Or(img.eq(40)) # 80=Water, 40=Veg
        vectors = img.updateMask(mask).reduceToVectors(geometry=aoi, scale=20, maxPixels=1e9)
        
        features = []
        for f in vectors.getInfo()['features']:
            lbl = f['properties']['label']
            f['properties']['custom_type'] = 'water' if lbl == 80 else 'veg'
            features.append(f)
            
        print(f"   -> GEE found {len(features)} items")
        return features
    except Exception as e:
        print(f"GEE Error: {e}"); return []

# --- FILE GENERATION ---
def generate_files(features, arp, radius_m):
    kml = etree.Element("kml", xmlns="http://www.opengis.net/kml/2.2")
    doc = etree.SubElement(kml, "Document")
    
    styles = {"water": "80ff0000", "veg": "80008000", "waste": "9913458b", "radius": "ff0000ff"}
    for k, c in styles.items():
        s = etree.SubElement(doc, "Style", id=k)
        etree.SubElement(etree.SubElement(s, "PolyStyle"), "color").text = c
        if k == 'radius': etree.SubElement(etree.SubElement(s, "LineStyle"), "color").text = c
    
    pm = etree.SubElement(doc, "Placemark")
    etree.SubElement(pm, "name").text = "ARP"
    etree.SubElement(etree.SubElement(pm, "Point"), "coordinates").text = f"{arp['lon']},{arp['lat']},0"
    
    csv_rows = []
    for f in features:
        t = f['properties'].get('custom_type', 'waste')
        area = f['properties'].get('area_sq_m', 0)
        name = "Water Body" if t == 'water' else "Vegetation" if t == 'veg' else "Industrial/Waste"
        
        pm = etree.SubElement(doc, "Placemark")
        etree.SubElement(pm, "name").text = name
        etree.SubElement(pm, "styleUrl").text = f"#{t}"
        etree.SubElement(pm, "description").text = f"Area: {int(area)}m2"
        
        geom = f['geometry']
        if geom['type'] in ['Polygon', 'MultiPolygon']:
            coords_list = geom['coordinates'][0] if geom['type'] == 'Polygon' else geom['coordinates'][0][0]
            coords_str = " ".join([f"{c[0]},{c[1]},0" for c in coords_list])
            poly = etree.SubElement(pm, "Polygon")
            etree.SubElement(etree.SubElement(poly, "outerBoundaryIs"), "LinearRing").append(etree.fromstring(f"<coordinates>{coords_str}</coordinates>"))

        c = shape(geom).centroid
        csv_rows.append({"Type": name, "Area_m2": int(area), "Distance_km": round(haversine(arp['lat'], arp['lon'], c.y, c.x), 2), "Lat": c.y, "Lon": c.x})

    return etree.tostring(kml).decode('utf-8'), pd.DataFrame(csv_rows).to_csv(index=False)

# --- ENDPOINT ---
@app.route('/generate-report', methods=['POST', 'OPTIONS'])
def generate_report():
    if request.method == 'OPTIONS': return jsonify({"status": "ok"}), 200
    try:
        d = request.json
        r_km = float(d.get('radius_km', 13))
        min_area = float(d.get('min_area_sq_m', 5000))
        
        if d.get('mode') == 'icao':
            q = f'[out:json];node["icao"="{d["icao"].upper()}"];out center;'
            res = requests.post("https://overpass-api.de/api/interpreter", data={'data': q}).json()['elements'][0]
            arp = {"name": d["icao"], "lat": res['lat'], "lon": res['lon']}
        else:
            arp = {"name": "Custom", "lat": float(d['lat']), "lon": float(d['lon'])}

        all_raw = fetch_osm_data(arp['lat'], arp['lon'], r_km*1000) + fetch_gee_data(arp['lat'], arp['lon'], r_km*1000)
        
        final, display = [], []
        for f in all_raw:
            s = shape(f['geometry'])
            area = calculate_area(f['geometry'])
            
            # Relaxation: If it's a point/line that we buffered, it's valid.
            # Otherwise apply min_area filter.
            if area < min_area and area > 100: continue 
            
            f['properties']['area_sq_m'] = area
            final.append(f)
            
            display.append({
                "type": "Feature", "properties": f['properties'],
                "geometry": mapping(s.simplify(0.001))
            })
            
        kml, csv = generate_files(final, arp, r_km*1000)
        return jsonify({
            "feature_count": len(final), "airport_info": arp,
            "map_geojson": {"type": "FeatureCollection", "features": display},
            "kml_string": kml, "csv_string": csv
        })
    except Exception as e:
        print(f"ERROR: {e}")
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5000)