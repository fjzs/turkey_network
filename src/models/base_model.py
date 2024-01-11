import abc

# https://www.scaler.com/topics/interface-in-python/

class BaseModel(abc.ABC):
    
    @abc.abstractmethod
    def solve(self):
        pass
    
    @abc.abstractmethod
    def save_solution(self):
        pass