from class_city import City
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Point
from typing import Dict


def plot_map(cities: Dict[int, City]):
    
    # Build df
    data = []
    for _, c in cities.items():
        data.append([c.longitude, c.latitude])
    df = pd.DataFrame(data, columns=['Longitude', 'Latitude'])
    
    # Print map
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude), crs="EPSG:4326") # turkey cities
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    world = world.query("name=='Turkey'")
    ax = world.plot(color="whitesmoke", edgecolor="black", linewidth=1.5)
    gdf.plot(ax=ax, color="red", markersize=10)
    plt.title("Cities of Turkey")
    plt.ylabel("Latitude")
    plt.xlabel("Longitude")

    # To plot connections between cities
    plt.plot([cities[1].longitude, cities[2].longitude], [cities[1].latitude, cities[2].latitude], color="black", linewidth=1)
    
    # To plot hubs
    plt.scatter(cities[1].longitude, cities[1].latitude, s=80, marker='s', color='b', alpha = 0.5)    

    plt.show()




    

if __name__ == "__main__":
    from dataloader import load_data
    cities_data = load_data()
    plot_map(cities_data)



