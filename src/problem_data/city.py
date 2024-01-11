from dataclasses import dataclass
from typing import Dict

@dataclass
class City:
    name: str
    id: int
    fixed_hub_cost: float
    latitude: float
    longitude: float
    
    # Data related to other cities, where this city is the origin
    distance_km_to_other_cities: Dict[int, float]
    flow_goods_to_other_cities: Dict[int, float]
    fixed_link_cost_to_other_cities: Dict[int, float]
    