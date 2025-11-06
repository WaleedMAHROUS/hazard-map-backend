import os
import time
import json
from geemap.conversion import ee_to_geojson
import ee
import requests
import geojson
import folium
from lxml import etree
import pandas as pd
from shapely.geometry import shape
from pyproj import Geod
from math import radians, cos, sin, asin, sqrt, degrees, atan2, pi
from flask import Flask, request, jsonify
from flask_cors import CORS

# -----------------------------------------------------------------
# STEP 1: GEE AUTHENTICATION & FLASK APP SETUP
# -----------------------------------------------------------------

SERVICE_ACCOUNT_KEY = 'service-account-key.json'

try:
    with open(SERVICE_ACCOUNT_KEY, 'r') as f:
        key_data = json.load(f)
        service_account_email = key_data.get('client_email')

    if not service_account_email:
        raise Exception("Could not find 'client_email' in the service-account-key.json file.")

    credentials = ee.ServiceAccountCredentials(
        service_account_email,
        SERVICE_ACCOUNT_KEY
    )
    ee.Initialize(credentials=credentials)
    print("GEE Authentication Successful.")

except Exception as e:
    print(f"GEE Authentication Failed: {e}")

app = Flask(__name__)
CORS(app)

# -----------------------------------------------------------------
# STEP 2: ALL HELPER FUNCTIONS
# -----------------------------------------------------------------

OVERPASS_API_URL = "https://overpass-api.de/api/interpreter"
WORKING_SERVER = None
HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/97.0.4692.71 Safari/537.36'
}

def get_bounding_box(lat, lon, radius_km):
    R = 6371
    lat_rad = radians(lat)
    lon_rad = radians(lon)
    lat_dist = radius_km / R
    min_lat = degrees(lat_rad - lat_dist)
    max_lat = degrees(lat_rad + lat_dist)
    lon_dist_rad = radius_km / (R * cos(lat_rad))
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

def get_feature_centroid(geom):
    s = shape(geom)
    centroid = s.centroid
    return centroid.y, centroid.x

def calculate_polygon_area(geom):
    try:
        s = shape(geom)
        if s.geom_type != 'Polygon':
            return 0
        geod = Geod(ellps='WGS84')
        area, _ = geod.geometry_area_perimeter(s)
        return abs(area)
    except Exception:
        return 0

def fetch_icao_data(icao_code):
    global WORKING_SERVER 
    WORKING_SERVER = None 
    print(f"Fetching data for ICAO: {icao_code}...")
    query = f"""
        [out:json][timeout:10];
        (
          node["icao"="{icao_code.upper()}"];
          way["icao"="{icao_code.upper()}"];
          relation["icao"="{icao_code.upper()}"];
        );
        out center;
    """
    server_url = OVERPASS_API_URL
    print(f"  -> Trying server: {server_url}")
    try:
        response = requests.post(server_url, data={'data': query}, headers=HEADERS, timeout=20)
        response.raise_for_status() 
        data = response.json()
        if not data['elements']:
            raise Exception(f"ICAO code '{icao_code}' not found.")
        
        airport = data['elements'][0]
        name = airport.get('tags', {}).get('name', f"{icao_code.upper()} Airport")
        if 'lat' in airport and 'lon' in airport:
            lat, lon = airport['lat'], airport['lon']
        elif 'center' in airport:
            lat, lon = airport['center']['lat'], airport['center']['lon']
        else:
            raise Exception("Could not find coordinates for this airport.")
            
        print(f"Found: {name} ({lat}, {lon})")
        WORKING_SERVER = server_url 
        return {"name": name, "lat": lat, "lon": lon} 
    except Exception as e:
        print(f"     ...An error occurred on this server: {e}")
        raise e

def fetch_hazard_data(lat, lon, radius_meters):
    global WORKING_SERVER 
    if not WORKING_SERVER:
        print("Warning: No working Overpass server. Skipping OSM.")
        return None
    
    radius_km_with_buffer = (radius_meters / 1000) + 5
    bbox = get_bounding_box(lat, lon, radius_km_with_buffer)
    bbox_str = f"{bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}"
    print(f"Fetching OSM hazard data using bounding box: {bbox_str}")
    all_features = []
    queries = {
        "waste_sites": (
            f'nwr["landuse"="landfill"]({bbox_str});'
            f'nwr["industrial"="waste_disposal"]({bbox_str});'
            f'nwr["amenity"="waste_disposal"]({bbox_str});'
            f'nwr["amenity"="slaughterhouse"]({bbox_str});'
            f'nwr["man_made"="wastewater_plant"]({bbox_str});'
        ),
        "water_bodies": (
            f'nwr["natural"="water"]({bbox_str});'
            f'nwr["natural"="wetland"]({bbox_str});'
            f'nwr["landuse"="reservoir"]({bbox_str});'
        ),
        "waterways": (f'nwr["waterway"~"river|canal|drain"]({bbox_str});'),
        "farms": (f'nwr["landuse"~"farmland|orchard"]({bbox_str});'),
        "forests": (f'nwr["landuse"="forest"]({bbox_str});' f'nwr["natural"="wood"]({bbox_str});')
    }
    
    for name, q_part in queries.items():
        print(f"  -> Fetching OSM {name}...")
        query = f"[out:json][timeout:60]; ({q_part}); out geom;"
        try:
            response = requests.post(WORKING_SERVER, data={'data': query}, headers=HEADERS, timeout=70)
            response.raise_for_status()
            if not response.text:
                continue
            data = response.json()
            if 'features' in data:
                all_features.extend(data['features'])
        except Exception as e:
            print(f"     ...Warning: OSM query for {name} failed: {e}")
            pass 
        time.sleep(2) 
            
    if not all_features:
        print("No mappable hazard features were found in OpenStreetMap for this area.")
        return None
        
    print(f"Total OSM features loaded: {len(all_features)}")
    return {"type": "FeatureCollection", "features": all_features}

def fetch_lulc_data(center_lat, center_lon, radius_meters):
    print(f"Fetching LULC data from Google Earth Engine...")
    try:
        aoi = ee.Geometry.Point(center_lon, center_lat).buffer(radius_meters)
        world_cover = ee.ImageCollection('ESA/WorldCover/v100').first()
        water = world_cover.eq(80).multiply(1)
        cropland = world_cover.eq(40).multiply(2)
        hazards_image = water.add(cropland).selfMask()
        vector_features = hazards_image.reduceToVectors(
            geometry=aoi,
            scale=30,
            maxPixels=1e10,
            geometryType='polygon'
        )
        if vector_features.size().getInfo() == 0:
            print(" -> No LULC features found in GEE for this area.")
            return None
        print(" -> GEE data found, converting to GeoJSON...")
        geojson_data = ee_to_geojson(vector_features) # Use direct import
        formatted_features = []
        for feature in geojson_data.get('features', []):
            props = feature.get('properties', {})
            label = props.get('label')
            new_tags = {}
            if label == 1:
                new_tags['natural'] = 'water'
                new_tags['name'] = 'Water Body (from Satellite)'
            elif label == 2:
                new_tags['landuse'] = 'farmland'
                new_tags['name'] = 'Cropland (from Satellite)'
            else:
                continue 
            feature['properties']['tags'] = new_tags
            formatted_features.append(feature)
        print(f" -> Processed {len(formatted_features)} LULC features from GEE.")
        return {"type": "FeatureCollection", "features": formatted_features}
    except Exception as e:
        print(f"   ...ERROR fetching GEE data: {e}")
        raise Exception(f"Google Earth Engine query failed: {e}")

def get_hazard_style(tags):
    tags = tags or {} 
    if tags.get('natural') in ('water', 'wetland') or tags.get('landuse') == 'reservoir' or 'waterway' in tags:
        return {'styleUrl': '#style_water'}
    if tags.get('landuse') == 'landfill' or tags.get('amenity') == 'waste_disposal' or tags.get('man_made') == 'wastewater_plant':
        return {'styleUrl': '#style_landfill'}
    if tags.get('amenity') == 'slaughterhouse':
        return {'styleUrl': '#style_other'}
    if tags.get('landuse') in ('farmland', 'orchard', 'forest') or tags.get('natural') == 'wood':
        return {'styleUrl': '#style_vegetation'}
    return {'styleUrl': '#style_default'}

def generate_csv_string(airport_info, hazard_data, min_area_sq_meters=0):
    if not hazard_data:
        print("No hazard data to generate CSV.")
        return None
    
    print(f"Generating CSV string (filtering < {min_area_sq_meters} sq.m)...")
    arp_lat, arp_lon = airport_info['lat'], airport_info['lon']
    report_data = []
    features_skipped = 0

    for feature in hazard_data.get('features', []): 
        tags = feature.get('properties', {}).get('tags', {})
        geom = feature.get('geometry')
        if not geom: continue
            
        type_str = tags.get('landuse') or tags.get('natural') or tags.get('amenity') or tags.get('waterway') or 'unknown'
        name_str = tags.get('name', f"Unnamed {type_str}")
        centroid_lat, centroid_lon = get_feature_centroid(geom)
        distance_km = haversine_distance(arp_lat, arp_lon, centroid_lat, centroid_lon)
        
        area_sq_m = 0
        if geom.get('type') == 'Polygon':
            area_sq_m = calculate_polygon_area(geom)
            if area_sq_m < min_area_sq_meters:
                features_skipped += 1
                continue 
        
        report_data.append({
            "name": name_str,
            "type": type_str,
            "latitude": centroid_lat,
            "longitude": centroid_lon,
            "distance_from_arp_km": distance_km,
            "area_sq_meters": area_sq_m,
            "osm_id": feature.get('id', 'N/A')
        })

    print(f"     ...Skipped {features_skipped} polygons.")
    df = pd.DataFrame(report_data)
    df = df.sort_values(by="distance_from_arp_km")
    return df.to_csv(index=False, encoding='utf-8-sig')

def generate_kml_string(airport_info, hazard_data, radius_meters, min_area_sq_meters=0):
    print(f"Generating KML string (filtering < {min_area_sq_meters} sq.m)...")
    kml = etree.Element("kml", xmlns="http://www.opengis.net/kml/2.2")
    doc = etree.SubElement(kml, "Document")
    doc.append(etree.Comment(f"Generated for {airport_info['name']}"))

    def create_style(id, color, fill_opacity=0.3):
        style = etree.SubElement(doc, "Style", id=id)
        line_style = etree.SubElement(style, "LineStyle")
        line_style.append(etree.XML(f"<color>ff{color[4:6]}{color[2:4]}{color[0:2]}</color>"))
        line_style.append(etree.XML("<width>2</width>"))
        poly_style = etree.SubElement(style, "PolyStyle")
        poly_style.append(etree.XML(f"<color>{hex(int(fill_opacity * 255))[2:].zfill(2)}{color[4:6]}{color[2:4]}{color[0:2]}</color>"))
    
    style_circle = etree.SubElement(doc, "Style", id="style_circle")
    etree.SubElement(style_circle, "LineStyle").extend([etree.XML("<color>ff0000ff</color>"), etree.XML("<width>2</width>")])
    etree.SubElement(style_circle, "PolyStyle").extend([etree.XML("<color>1aff0000</color>"), etree.XML("<outline>1</outline>"), etree.XML("<fill>1</fill>")])
    create_style("style_water", "0000FF", 0.4)
    create_style("style_landfill", "8B4513", 0.5)
    create_style("style_vegetation", "228B22", 0.3)
    create_style("style_other", "800080", 0.5)
    create_style("style_default", "FFFF00", 0.3)
    
    pm = etree.SubElement(doc, "Placemark")
    etree.SubElement(pm, "name").text = f"{airport_info['name']} (ARP / Center)"
    etree.SubElement(etree.SubElement(pm, "Point"), "coordinates").text = f"{airport_info['lon']},{airport_info['lat']},0"
    
    pm_circle = etree.SubElement(doc, "Placemark")
    etree.SubElement(pm_circle, "name").text = f"{radius_meters/1000}km Radius"
    etree.SubElement(pm_circle, "styleUrl").text = "#style_circle"
    lat, lon = airport_info['lat'], airport_info['lon']
    R = 6378137
    d = radius_meters / R
    lat_rad, lon_rad = radians(lat), radians(lon)
    points = []
    for i in range(73):
        brng = radians(i * 5)
        new_lat_rad = asin(sin(lat_rad) * cos(d) + cos(lat_rad) * sin(d) * cos(brng))
        new_lon_rad = lon_rad + atan2(sin(brng) * sin(d) * cos(lat_rad), cos(d) - sin(lat_rad) * sin(new_lat_rad))
        points.append(f"{degrees(new_lon_rad)},{degrees(new_lat_rad)},0")
    etree.SubElement(etree.SubElement(etree.SubElement(pm_circle, "Polygon"), "outerBoundaryIs"), "LinearRing").append(etree.XML(f"<coordinates>{' '.join(points)}</coordinates>"))

    features_skipped = 0
    features_kept = 0
    if hazard_data:
        for feature in hazard_data.get('features', []):
            try:
                props = feature.get('properties', {})
                tags = props.get('tags', {})
                name = tags.get('name', props.get('type', 'Hazard'))
                style_url = get_hazard_style(tags)['styleUrl']
                geom = feature['geometry']
                if not geom: continue

                area_sq_m = calculate_polygon_area(geom)
                if geom.get('type') == 'Polygon' and area_sq_m < min_area_sq_meters:
                    features_skipped += 1
                    continue
                
                features_kept += 1
                pm_haz = etree.SubElement(doc, "Placemark")
                etree.SubElement(pm_haz, "name").text = name
                etree.SubElement(pm_haz, "styleUrl").text = style_url
                
                centroid_lat, centroid_lon = get_feature_centroid(geom)
                distance_km = haversine_distance(airport_info['lat'], airport_info['lon'], centroid_lat, centroid_lon)
                desc = f"<b>Distance from ARP:</b> {distance_km:.2f} km<br>"
                if area_sq_m > 0: desc += f"<b>Area:</b> {area_sq_m:,.0f} sq. meters<br><hr>"
                for k, v in tags.items(): desc += f"<b>{k}</b>: {v}<br>"
                etree.SubElement(pm_haz, "description").text = etree.CDATA(desc)
            
                if geom['type'] == 'Point':
                    etree.SubElement(etree.SubElement(pm_haz, "Point"), "coordinates").text = f"{geom['coordinates'][0]},{geom['coordinates'][1]},0"
                elif geom['type'] == 'Polygon':
                    coords = " ".join([f"{c[0]},{c[1]},0" for c in geom['coordinates'][0]])
                    etree.SubElement(etree.SubElement(etree.SubElement(pm_haz, "Polygon"), "outerBoundaryIs"), "LinearRing").append(etree.XML(f"<coordinates>{coords}</coordinates>"))
                elif geom['type'] == 'LineString':
                    coords = " ".join([f"{c[0]},{c[1]},0" for c in geom['coordinates']])
                    etree.SubElement(etree.SubElement(pm_haz, "LineString"), "coordinates").text = coords
            except Exception:
                continue
    
    print(f"     ...Skipped {features_skipped} polygons.")
    print(f"     ...Added {features_kept} features to KML.")
    return etree.tostring(kml, pretty_print=True, xml_declaration=True, encoding='UTF-8')


# --- REMOVED STEP 3 (COMPASS ANALYSIS) ---


# -----------------------------------------------------------------
# STEP 4: THE API ENDPOINT (Simplified)
# -----------------------------------------------------------------

@app.route('/generate-report', methods=['POST', 'OPTIONS'])
def generate_report():
    if request.method == 'OPTIONS':
        return jsonify({"status": "ok"}), 200

    print("\n--- Received new request ---")
    
    try:
        data = request.json
        mode = data.get('mode')
        radius_km = float(data.get('radius_km', 13))
        min_area_sq_m = float(data.get('min_area_sq_m', 0))
        HAZARD_RADIUS_METERS = radius_km * 1000

        airport_info = None

        if mode == 'icao':
            icao = data.get('icao')
            if not icao:
                raise Exception("Mode is 'icao' but no 'icao' code was provided.")
            airport_info = fetch_icao_data(icao)
        elif mode == 'coords':
            lat = float(data.get('lat'))
            lon = float(data.get('lon'))
            if lat is None or lon is None:
                raise Exception("Mode is 'coords' but 'lat' or 'lon' was not provided.")
            airport_info = {"name": "Custom Location", "lat": lat, "lon": lon}
        else:
            raise Exception(f"Invalid 'mode' provided: {mode}")

        all_hazard_features = []
        
        osm_data = fetch_hazard_data(airport_info['lat'], airport_info['lon'], HAZARD_RADIUS_METERS)
        if osm_data and 'features' in osm_data:
            all_hazard_features.extend(osm_data['features'])

        lulc_data = fetch_lulc_data(airport_info['lat'], airport_info['lon'], HAZARD_RADIUS_METERS)
        if lulc_data and 'features' in lulc_data:
            all_hazard_features.extend(lulc_data['features'])

        hazard_data = None
        if all_hazard_features:
            hazard_data = {"type": "FeatureCollection", "features": all_hazard_features}
            print(f"Total features found (before filtering): {len(all_hazard_features)}")
        else:
            print("No hazards found in OSM or GEE LULC data.")

        kml_string = generate_kml_string(airport_info, hazard_data, HAZARD_RADIUS_METERS, min_area_sq_m)
        csv_string = generate_csv_string(airport_info, hazard_data, min_area_sq_m)

        features_for_map = []
        if hazard_data:
            for f in hazard_data['features']:
                geom = f.get('geometry')
                if not geom: continue
                area_sq_m = 0 
                if geom.get('type') == 'Polygon':
                    area_sq_m = calculate_polygon_area(geom)
                if geom.get('type') != 'Polygon' or area_sq_m >= min_area_sq_m:
                    features_for_map.append(f)
        
        map_geojson = {"type": "FeatureCollection", "features": features_for_map}
        feature_count = len(features_for_map)
        
        # --- ANALYSIS REMOVED ---

        print(f"--- Request complete. Returning {feature_count} features. ---")

        # 6. Send all data back to the frontend
        return jsonify({
            "message": f"Success! Found {feature_count} features.", # Simpler message
            "feature_count": feature_count,
            "airport_info": airport_info,
            "kml_string": kml_string.decode('utf-8'),
            "csv_string": csv_string,
            "map_geojson": map_geojson,
            "analysis_summary": None # Send null
        }), 200

    except Exception as e:
        print(f"--- Request FAILED: {e} ---")
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5000)