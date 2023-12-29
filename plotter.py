import geopandas as gpd
import matplotlib.pyplot as plt
import os
from utils import FOLDER


def plot_map(cities_data: dict(), 
             cities_considered: list, 
             hubs_ids: set(), 
             collection: dict(), 
             transfer: dict(), 
             distribution: dict(),
             plot_name: str):
    """Plots the map and optionally the solution

    Args:
        cities (Dict[int, City]): data of the cities
        cities_considered (list): ids of cities to plot
        solution_hubs_ids (set of int): set of hubs
        solution_city_hub_connection (set of int, int): set of (i, j) where i is a city and j is a hub

    Args:
        cities (dict): data of the cities
        hubs_ids (set): ids of the hubs
        collection (dict): collection flow: (i,k) -> flow
        transfer (dict): transfer flow: (i,k,l) -> flow
        distribution (dict): distribution flow (i,l,j) -> flow
    """
    
    # Print map
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    world = world.query("name == 'Turkey'")
    ax = world.plot(color="whitesmoke", edgecolor="black", linewidth=1.5)
    ax.figure.set_size_inches(14,6)
    plt.title(plot_name)
    plt.ylabel("Latitude")
    plt.xlabel("Longitude")

    # Plot cities id, but if a solution is given, only plot those in the solution
    cities_to_plot = {id for id in cities_data.keys()}
    if cities_considered is not None:
        cities_to_plot = cities_considered
    
    # Plot the cities numbers
    for i in cities_to_plot:
        city = cities_data[i]
        plt.annotate(str(i), (city.longitude-0.1, city.latitude-0.1), color="white", size=10, zorder = 100)
        
        # Plot hub shape or regular city
        if hubs_ids is not None:
            if i in hubs_ids:
                plt.scatter(cities_data[i].longitude, 
                            cities_data[i].latitude, 
                            s=120, 
                            marker='s', 
                            color='blue', 
                            zorder=80)
            else:
                plt.scatter(cities_data[i].longitude, 
                            cities_data[i].latitude, 
                            s=120, 
                            marker="o", 
                            color="red", 
                            zorder=80)

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

    # Get max flow to gauge linewidths
    max_flow = max(total_flow.values())
    min_flow = min(total_flow.values())
    min_linewidth = 0.5
    max_linewidth = 8
    def get_linewidth_given_flow(x) -> float:
        return  min_linewidth + ((x - min_flow)/(max_flow - min_flow + 1)) * (max_linewidth - min_linewidth)

    # Plot flows between cities
    for (i,j), flow in total_flow.items():
        plt.plot([cities_data[i].longitude, cities_data[j].longitude], 
                 [cities_data[i].latitude,   cities_data[j].latitude], 
                 color="black", 
                 linewidth=get_linewidth_given_flow(flow),
                 zorder=10)
    
    # Save plot
    filepath = os.path.join(FOLDER, plot_name + ".png")
    plt.savefig(filepath)
    print(f"\nPlot saved in {filepath}")
    

# if __name__ == "__main__":
#     from dataloader import load_data
#     cities_data = load_data()
#     plot_map(cities_data, solution_hubs_ids=None, solution_city_hub_connection=None)



