from .city import City
from .dataloader import load_data
from typing import Dict

class ProblemData:
    """This class is responsible for providing the data needed for
    this problem that is located in files. In this case, its just
    a dictionary of cities, in the form: {id: City}
    """
    
    def __init__(self):
        self.cities = load_data()
        