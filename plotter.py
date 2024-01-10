from class_city import City
import geopandas as gpd
import matplotlib.pyplot as plt
import os
from typing import Dict
from utils import FOLDER

MIN_LINEWIDTH = 0.5
MAX_LINEWIDTH = 4
MAP_LINEWIDTH = 1.5
MIN_SIZEBALL = 5
MAX_SIZEBALL = 600

def _get_size_given_min_max(x, min_x, max_x, min_size, max_size) -> float:
    """This just creates the uniform distribution

    Args:
        x:
        min_x:
        max_x:
        min_size:
        max_size:

    Returns:
        float:
    """
    return  min_size + ((x - min_x)/(max_x - min_x + 0.01)) * (max_size - min_size)

def plot_map(cities_data: Dict[int, City],
             cities_considered: list, 
             hubs_ids: set(), 
             collection: dict(), 
             transfer: dict(), 
             distribution: dict(),
             beta: float,
             alpha: float,
             top_k_cities_for_hub: set,
             size_proportional_to_flow: bool,
             draw_city_names: bool):
    """Plots the map and optionally the solution

    Args:
        - cities_data (dict): data of the cities
        - cities_considered (list): cities to plot
        - hubs_ids (set): ids of the hubs
        - collection (dict): collection flow: (i,k) -> flow
        - transfer (dict): transfer flow: (i,k,l) -> flow
        - distribution (dict): distribution flow (i,l,j) -> flow
        - beta (float): max delivery time between origin and destination
        - alpha (float): hub-to-hub cost factor (greater than truck transportation)
        - top_k_cities_for_hub (list): top k cities to be assigned as potential hubs
        - size_proportional_to_flow (bool): to see the map and the volumes required to transport
        - draw_city_names (bool)
    """
    
    assert isinstance(top_k_cities_for_hub, set)

    # Print map
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    world = world.query("name == 'Turkey'")
    ax = world.plot(color="whitesmoke", edgecolor="black", linewidth=MAP_LINEWIDTH)
    ax.figure.set_size_inches(14,6)

    plt.title(f"N={len(cities_considered)}, β={beta}[h], α={alpha}, k={top_k_cities_for_hub}")
    plt.ylabel("Latitude")
    plt.xlabel("Longitude")

    # Plot cities id, but if a solution is given, only plot those in the solution
    cities_to_plot = {id for id in cities_data.keys()}
    if cities_considered is not None:
        cities_to_plot = cities_considered
    
    # Check volumes
    if size_proportional_to_flow:
        supply_per_city = dict()
        for i, city in cities_data.items():
            supply_per_city[i] = sum(city.flow_goods_to_other_cities.values())
        min_supply = min(supply_per_city.values())
        max_supply = max(supply_per_city.values())

    # Plot the cities numbers
    for i in cities_to_plot:
        city = cities_data[i]
        
        # Plotting a solution
        # if hubs_ids is not None:
        #     plt.annotate(str(i), (city.longitude-0.1, city.latitude-0.1), color="white", size=10, zorder=100)
        #     if i in hubs_ids:
        #         plt.scatter(cities_data[i].longitude, cities_data[i].latitude, s=150, marker='s', color='blue', zorder=80)
        #     else:
        #         node_color = "red" if i in non_hubs_nodes else "green"
        #         plt.scatter(cities_data[i].longitude + 0.05, cities_data[i].latitude, s=160, marker="o", color=node_color, zorder=80)
        
        # Plotting proportional to volume per node
        
        size_ball = 100
        if size_proportional_to_flow:
            size_ball = _get_size_given_min_max(supply_per_city[i], min_supply, max_supply, min_size=MIN_SIZEBALL, max_size=MAX_SIZEBALL)
        
        node_color = "blue" if i in hubs_ids else ("red" if i in top_k_cities_for_hub else "green")
        plt.scatter(cities_data[i].longitude, 
                    cities_data[i].latitude, 
                    s=size_ball, 
                    marker="o", 
                    color=node_color, 
                    zorder=80)
        
        if draw_city_names:
            plt.annotate(str(city.name)[0:3], (city.longitude, city.latitude), color="black", size=10, zorder=100)
        
        # Plotting just the numbers
        #plt.scatter(cities_data[i].longitude + 0.05, cities_data[i].latitude, s=160, marker="o", color="red", zorder=80)

    # Compute the total flow between each pair of cities (sum (i,j) and (j,i) in the same arc)
    total_flow = dict() # (i,j) -> flow
    
    # Aggregate collection
    if collection is not None:
        for (i,k), flow in collection.items():
            if (i,k) not in total_flow:
                total_flow[(i,k)] = 0
            if (k,i) not in total_flow:
                total_flow[(k,i)] = 0
            total_flow[(i,k)] += flow
            total_flow[(k,i)] += flow

    # Aggregate transfers
    if transfer is not None:
        for (i,k,l), flow in transfer.items():
            if (k,l) not in total_flow:
                total_flow[(k,l)] = 0
            if (l,k) not in total_flow:
                total_flow[(l,k)] = 0
            total_flow[(k,l)] += flow
            total_flow[(l,k)] += flow
    
    # Aggregate distribution
    if distribution is not None:
        for (i,l,j), flow in distribution.items():
            if (l,j) not in total_flow:
                total_flow[(l,j)] = 0
            if (j,l) not in total_flow:
                total_flow[(j,l)] = 0
            total_flow[(l,j)] += flow
            total_flow[(j,l)] += flow

    # Remove redundant keys
    total_flow = {(i,j):v for (i,j),v in total_flow.items() if i < j}
    print(f"\nTotal flow in arcs:")
    for (i,j), flow in total_flow.items():
        print(f"{i,j}: {flow}")

    # Get max flow to gauge linewidths
    max_flow = 0 if len(total_flow) == 0 else max(total_flow.values())
    min_flow = 0 if len(total_flow) == 0 else min(total_flow.values())    

    # Plot flows between cities
    for (i,j), flow in total_flow.items():
        plt.plot([cities_data[i].longitude, cities_data[j].longitude], 
                 [cities_data[i].latitude,   cities_data[j].latitude], 
                 color="black", 
                 linewidth=  _get_size_given_min_max(flow, min_flow, max_flow, MIN_LINEWIDTH, MAX_LINEWIDTH),
                 zorder=10)
    
    # Save plot
    plot_name = "n" + str(len(cities_considered)).zfill(2) 
    plot_name += "_t" + str(beta).zfill(2) 
    plot_name += "_a" + str(alpha).zfill(2) 
    plot_name += "_h" + str(len(top_k_cities_for_hub)).zfill(2)
    
    filepath = os.path.join(FOLDER, plot_name + ".png")
    plt.savefig(filepath)
    print(f"\nPlot saved in {filepath}")
    

if __name__ == "__main__":
    from dataloader import load_data
    cities_data = load_data()
    plot_map(cities_data,
             cities_considered=None,
             hubs_ids=None,
             collection=None,
             transfer=None,
             distribution=None,
             plot_name="test",
             size_proportional_to_flow=True)



