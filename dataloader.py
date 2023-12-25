from class_city import City
import pandas as pd
from typing import Dict

DATA_FILE = "./data/data.xlsx"


def __get_cities(data_file:str = DATA_FILE):
    """Returns an initial dictionary of cities with information pertaining
    only to the city

    Args:
        data_file (str, optional): Defaults to DATA_FILE.

    Returns:
        cities (dictionary)
    """

    df = pd.read_excel(data_file, sheet_name="cities", header=0, index_col=0)
    cities = dict()

    for row in df.iterrows():
        id = row[0]
        name = row[1]["City"]
        latitude = row[1]["Latitude"]
        longitude = row[1]["Longitude"]
        fixed_hub_cost = row[1]["Fixed Hub Cost"]
        city = City(name, id, fixed_hub_cost, latitude, longitude, dict(), dict(), dict(), dict())
        cities[id] = city
    
    return cities


def __fill_city_to_city_info(cities: Dict[int, City], sheet_name: str, data_file:str = DATA_FILE):
    """Fills one parameter related to a city-to-city relationship

    Args:
        cities (Dict[int, City]): 
        sheet_name (str): 
        data_file (str, optional):Defaults to DATA_FILE.

    Raises:
        ValueError:
    
    Returns:
        None
    """

    df = pd.read_excel(data_file, sheet_name=sheet_name, header=1, index_col=1)
    for row in df.iterrows():
        origin_id = row[0]
        info_to_other_city = row[1].to_dict() # {'Name': 'ADANA', 1: 0, 2: 329, 3: 573, 4: 961, ...
        for dest_id, value in info_to_other_city.items():
            if dest_id != "Name":
                if sheet_name == "distance_km":
                    cities[origin_id].distance_km_to_other_cities[dest_id] = value
                elif sheet_name == "travel_time_min":
                    cities[origin_id].travel_time_min_to_other_cities[dest_id] = value
                elif sheet_name == "flow":
                    cities[origin_id].flow_goods_to_other_cities[dest_id] = value
                elif sheet_name == "fixed_link_cost":
                    cities[origin_id].fixed_link_cost_to_other_cities[dest_id] = value
                else:
                    raise ValueError(f"I dont know this sheet name: {sheet_name}")
        


def load_data() -> Dict[int, City]:
    """Loads the data for this problem

    Returns:
        cities (Dict[int, City]): dictionary that maps every city id to its object
    """
    cities = __get_cities()
    __fill_city_to_city_info(cities, sheet_name = "distance_km")
    __fill_city_to_city_info(cities, sheet_name = "travel_time_min")
    __fill_city_to_city_info(cities, sheet_name = "flow")
    __fill_city_to_city_info(cities, sheet_name = "fixed_link_cost")
    return cities


if __name__ == "__main__":
    load_data()
    
