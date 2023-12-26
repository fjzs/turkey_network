from class_city import City
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from typing import Dict


def plot_map(cities: Dict[int, City], solution_hubs_ids: set(), solution_city_hub_connection: set()):
    """Plots the map and optionally the solution

    Args:
        cities (Dict[int, City]): data of the cities
        solution_hubs_ids (set of int): set of hubs
        solution_city_hub_connection (set of int, int): set of (i, j) where i is a city and j is a hub
    """
    
    # Build df
    data = []
    for _, c in cities.items():
        data.append([c.longitude, c.latitude])
    df = pd.DataFrame(data, columns=['Longitude', 'Latitude'])
    
    # Print map
    #gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude), crs="EPSG:4326")
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    world = world.query("name == 'Turkey'")
    ax = world.plot(color="whitesmoke", edgecolor="black", linewidth=1.5)
    #gdf.plot(ax=ax, color="red", markersize=10)
    ax.figure.set_size_inches(14,6)
    plt.title("Cities of Turkey")
    plt.ylabel("Latitude")
    plt.xlabel("Longitude")

    # To plot cities
    for id, city in cities.items():
        plt.annotate(str(id), (city.longitude-0.1, city.latitude-0.1), color="red", size=10)

    # To plot hubs
    if solution_hubs_ids is not None:
        for i in solution_hubs_ids:
            plt.scatter(cities[i].longitude, cities[i].latitude, s=100, marker='s', color='grey', alpha = 0.7)
            #plt.annotate(str(i), (cities[i].longitude-0.1, cities[i].latitude-0.1))

    # To plot connections between cities
    if solution_city_hub_connection is not None:
        for (i,j) in solution_city_hub_connection:
            plt.plot([cities[i].longitude, cities[j].longitude], [cities[i].latitude, cities[j].latitude], color="black", linewidth=1)
    
    
    plt.show()
    

# if __name__ == "__main__":
#     from dataloader import load_data
#     cities_data = load_data()
#     plot_map(cities_data, solution_hubs_ids=None, solution_city_hub_connection=None)



