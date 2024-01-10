from dataclasses import dataclass

@dataclass
class City:
    name: str
    id: int
    fixed_hub_cost: float
    latitude: float
    longitude: float
    
    # Data related to other cities, where this city is the origin
    distance_km_to_other_cities: dict()
    travel_time_min_to_other_cities: dict()
    flow_goods_to_other_cities: dict()
    fixed_link_cost_to_other_cities: dict()
    