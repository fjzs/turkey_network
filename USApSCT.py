from typing import Dict
from class_city import City
import gurobipy as gp
from gurobipy import GRB
import plotter
import utils

class USAPSCT:
    """This class implements the uncapacitated, single allocation p-hub set covering with
    time service constraints, a mixture between:
    - Ernst, A.T., Krishnamoorthy, M. (1996)
    - Tan, P.Z., Kara, B. Y. (2007)
    
    References:
    - Tan, P.Z., Kara, B. Y. (2007). “A hub covering model for cargo delivery systems” 
    Networks 49(1): 28-39, doi: 10.1002/net.20139.
    - Ernst, A.T., Krishnamoorthy, M. (1996). “Efficient algorithms for the uncapacitated 
    ingle allocation p-hub median problem”. Location Science 4(3): 139-154.
    """

    def __init__(self, cities_data: Dict[int, City], 
                 max_nodes:int = 81, 
                 max_arrival_time_h:int = 30,
                 economy_of_scale_factor:float = 0.9,
                 top_k_cities_for_hub:int = 20,
                 min_latitude_to_be_hub: float = 36,
                 max_latitude_to_be_hub: float = 45,
                 min_longitude_to_be_hub: float = 26,
                 max_longitude_to_be_hub: float = 45
                 ):
        assert 2 <= max_nodes <= 81
        assert 1 <= max_arrival_time_h
        self.cities_data = cities_data
        
        # TODO:
        # Done: Eliminate infeasible combinations de i -> k -> l -> j because travel time is violated
        # Not needed: Optimize second phase objective to minimize travel times
        # Check cuts with presolve
        # https://github.com/csymvoul/python-structure-template

        # Define sets starting with set_
        self.set_cities = self._get_set_cities(max_nodes)                                   # N

        # Define parameters starting with par_
        self.par_flow_box = self._get_city_to_city_parameter("flow")                        # w_ij
        self.par_fixed_hub_cost = self._get_hub_cost()                                      # h_i
        self.par_distance_km = self._get_city_to_city_parameter("distance_km")              # d_ij
        self.par_fixed_link_cost = self._get_city_to_city_parameter("fixed_link_cost")      # f_ij
        self.par_supply = self._get_node_supply()                                           # s_i
        self.par_unit_cost_per_km_flow = 1                                                  # c
        self.par_speed_trucks_kph = 60                                                      # v
        self.par_travel_time_truck_h = self._get_travel_time(self.par_speed_trucks_kph)     # t_ij
        self.par_economy_scale_factor = economy_of_scale_factor                             # α
        self.par_max_arrival_time_h = max_arrival_time_h                                    # β
        
        # Heuristic inputs
        self.top_k_cities_for_hub = top_k_cities_for_hub
        self.min_latitude_to_be_hub = min_latitude_to_be_hub
        self.max_latitude_to_be_hub = max_latitude_to_be_hub
        self.min_longitude_to_be_hub = min_longitude_to_be_hub
        self.max_longitude_to_be_hub = max_longitude_to_be_hub
        self.heuristic_non_hubs_nodes = self._get_heuristic_non_hubs_nodes()

        # Create optimization model
        self.instance_name = "n" + str(len(self.set_cities)).zfill(2)
        self.instance_name += "_t" + str(self.par_max_arrival_time_h).zfill(2)
        self.instance_name += "_a" + str(self.par_economy_scale_factor).zfill(2)
        self.instance_name += "_k" + str(self.top_k_cities_for_hub).zfill(2) 
        print(f"Instance name: {self.instance_name}")
        self.model = gp.Model('USApSCT')
        
        # Create variables (var_)
        self.var_flow_orig_hub = self._add_var_flow_orig_hub()                              # F_ik
        self.var_assign_orig_hub = self._add_var_assign_orig_hub()                          # Z_ik
        self.var_flow_orig_hub_hub = self._add_var_flow_orig_hub_hub()                      # Y_ikl
        self.var_flow_orig_hub_dest = self._add_var_flow_orig_hub_dest()                    # X_ilj
        self.var_link_hubs = self._add_var_link_hubs()                                      # R_kl
        self.var_departure_from_hub_to_hub = self._add_var_departure_from_hub_to_hub()      # Dh_k
        self.var_departure_from_hub_to_dest = self._add_var_departure_from_hub_to_dest()    # Dc_k
        self.var_arrival_to_dest = self._add_var_arrival_to_dest()                          # A_j

        # Create constraints
        self._add_constraint_assign_node_to_single_hub()
        self._add_constraint_hub_implication()
        self._add_constraint_flow_from_source()
        self._add_constraint_flow_between_hubs_implies_link()
        self._add_constraint_conservation_of_flow_at_hub_from_origin()
        self._add_constraint_flow_node_in_hubs_implies_node_to_hub()
        self._add_constraint_flow_origin_destination_covered()
        self._add_constraint_link_between_hubs_implies_hubs()
        self._add_constraint_departure_from_hub_to_hub_waits_customers()
        self._add_constraint_departure_from_hub_to_destination_waits_for_hub_trucks()
        self._add_constraint_departure_from_hub_to_destination_waits_for_customer_trucks()
        self._add_constraint_latest_arrival_time_to_destination()
        self._add_constraint_latest_arrival_to_dest_within_limit()

        # Heuristics for solving big instances
        self._add_contraint_heuristic_non_hubs()

        # Set objective function
        # Variables for objective function (var_OF_)
        self.var_OF_location_cost = self._add_var_OF_location_cost()
        self.var_OF_transport_collection_cost = self._add_var_OF_transport_collection_cost()
        self.var_OF_transport_transfer_cost = self._add_var_OF_transport_transfer_cost()
        self.var_OF_transport_distribution_cost = self._add_var_OF_transport_distribution_cost()
        self.var_OF_links_cost = self._add_var_OF_links_cost()
        # Constraints for setting the variables of the OF
        self._add_constraint_OF_location_cost()
        self._add_constraint_OF_transport_collection_cost()
        self._add_constraint_OF_transport_transfer_cost()
        self._add_constraint_OF_transport_distribution_cost()
        self._add_constraint_OF_links_cost()
        # Setting the objecting function in the model
        self._set_objective_function()

        # Define a start value for some variables to speed up the process
        #self._set_start_values()

        self.model.update()

        # print(f"\n\nCONSTRAINTS:")
        # for i, con in enumerate(self.model.getConstrs()):
        #     print(f"\n{con}:\n{self.model.getRow(con)} {con.Sense} {con.RHS}")            
        

    
    def _get_set_cities(self, max_nodes: int) -> set:
        """Returns the set of cities to work with given the max_nodes input

        Args:
            max_nodes (int): max number of nodes to work with.

        Returns:
            set:
        """
        return set(list(cities_data.keys())[0: max_nodes])

    def _get_heuristic_non_hubs_nodes(self) -> set:
        """We limit here the possibility of some nodes to be hubs by
        - size of the city
        - min/max longitude
        - min/max latitude

        Returns:
            list: the non hubs nodes
        """
        non_hubs_nodes = set()
        print("\nApplying heuristic:")
        # The hubs with ranking > min_city_supply_ranking_to_be_hub wont be hubs 
        if self.top_k_cities_for_hub is not None:
            supply_ids = [(sum(city.flow_goods_to_other_cities.values()), id) for id, city in self.cities_data.items()]
            supply_ids.sort(reverse=True) # sorted in descending order
            print(f"Ordered cities in total supply:")
            for i, (val,id) in enumerate(supply_ids):
                name = self.cities_data[id].name
                print(f"# {i+1} has id {id}: {name} has total supply of {val}")

            # Remove every node which is not in the Top K
            hubs_to_remove = supply_ids[self.top_k_cities_for_hub:]
            hubs_to_remove = [id for (supply,id) in hubs_to_remove if id in self.set_cities]
            non_hubs_nodes.update(hubs_to_remove)

        # Remove eastern nodes given by min longitude
        if self.min_longitude_to_be_hub is not None:
            hubs_to_remove = [id for id in self.set_cities 
                              if self.cities_data[id].longitude < self.min_longitude_to_be_hub]
            non_hubs_nodes.update(hubs_to_remove)
        
        # Remove western nodes given by max longitude
        if self.max_longitude_to_be_hub is not None:
            hubs_to_remove = [id for id in self.set_cities 
                              if self.cities_data[id].longitude > self.max_longitude_to_be_hub]
            non_hubs_nodes.update(hubs_to_remove)
        
        # Remove northern nodes given by max latitude
        if self.max_latitude_to_be_hub is not None:
            hubs_to_remove = [id for id in self.set_cities 
                              if self.cities_data[id].latitude > self.max_latitude_to_be_hub]
            non_hubs_nodes.update(hubs_to_remove)
        
        # Remove southern nodes given by min latitude
        if self.min_latitude_to_be_hub is not None:
            hubs_to_remove = [id for id in self.set_cities 
                              if self.cities_data[id].latitude < self.min_latitude_to_be_hub]
            non_hubs_nodes.update(hubs_to_remove)
        
        print(f"\nHeuristic will remove {len(non_hubs_nodes)} nodes as potential hubs: {non_hubs_nodes}")

        return non_hubs_nodes


    def _get_city_to_city_parameter(self, param_name: str) -> dict():
        param = dict() # (i,j) -> parameter_ij
        for i in self.set_cities:
            city_i = self.cities_data[i]
            
            city_i_to_j_dict = dict()
            if param_name == "distance_km":
                city_i_to_j_dict = city_i.distance_km_to_other_cities
            elif param_name == "travel_time_min":
                city_i_to_j_dict = city_i.travel_time_min_to_other_cities
            elif param_name == "flow":
                city_i_to_j_dict = city_i.flow_goods_to_other_cities
            elif param_name == "fixed_link_cost":
                city_i_to_j_dict = city_i.fixed_link_cost_to_other_cities
            else:
                raise ValueError(f"param_name {param_name} not recognized")

            for j, value_ij in city_i_to_j_dict.items():
                if j in self.set_cities:
                    param[(i, j)] = value_ij
        return param

    def _get_hub_cost(self) -> dict():
        return {i: city.fixed_hub_cost for i, city in self.cities_data.items()}

    def _get_node_supply(self) -> dict:
        flow_origin = dict() # supply of node i
        for i in self.set_cities:
            city_i = self.cities_data[i]
            total_supply = 0
            for j, value_ij in city_i.flow_goods_to_other_cities.items():
                if j in self.set_cities:
                    total_supply += value_ij
            flow_origin[i] = total_supply
        return flow_origin

    def _get_travel_time(self, vehicle_speed: float) -> dict:
        assert vehicle_speed > 0
        travel_time = dict() # (i,j) -> parameter_ij
        for i in self.set_cities:
            city_i = self.cities_data[i]
            for j, dist_km_ij in city_i.distance_km_to_other_cities.items():
                travel_time[(i,j)] = dist_km_ij / vehicle_speed
        return travel_time


    #########################################
    #                                       #
    #              VARIABLES                #
    #                                       #
    #########################################

    def _add_var_flow_orig_hub(self):
        """
        F_ik: Flow from node i to hub in node k (>= 0)
        """
        N = self.set_cities
        return self.model.addVars(N, N, lb=0, vtype = GRB.CONTINUOUS, name="F")

    def _add_var_assign_orig_hub(self):
        """
        Z_ik: 1 if node i is assigned to hub at node k, else 0 (binary)
        """
        N = self.set_cities
        return self.model.addVars(N, N, lb=0, ub=1, vtype = GRB.BINARY, name="Z")

    def _add_var_flow_orig_hub_hub(self):
        """
        * Y_ikl: flow from i that goes between hubs k and l
        * Ɐ i ∈ N, Ɐ k ∈ N, Ɐ l ∈ N: k != l & i != l
        """
        N = self.set_cities
        N3 = set()
        for i in N:
            for k in N:
                for l in N:
                    if (k != l) and (i != l):
                        N3.add((i,k,l))
        return self.model.addVars(N3, lb=0, vtype = GRB.CONTINUOUS, name="Y")
    
    def _add_var_flow_orig_hub_dest(self):
        """
        * X_ilj: flow from i going through hub l to destination j
        * Ɐ i ∈ N, l ∈ N, j ∈ N: i != j
        """
        N = self.set_cities
        N3 = set()
        for i in N:
            for l in N:
                for j in N:
                    if i != j:
                        N3.add((i,l,j))
        return self.model.addVars(N3, lb=0, vtype = GRB.CONTINUOUS, name="X")

    def _add_var_link_hubs(self):
        """R_kl: 1 if hubs k and l are connected to transport goods
        """
        N = self.set_cities
        return self.model.addVars(N, N, lb=0, vtype = GRB.BINARY, name="R")

    def _add_var_departure_from_hub_to_hub(self):
        """D^h_k: departure time from hub k for transportation to other hubs
        """
        N = self.set_cities
        return self.model.addVars(N, lb=0, ub=self.par_max_arrival_time_h, vtype = GRB.CONTINUOUS, name="Dh")

    def _add_var_departure_from_hub_to_dest(self):
        """D^c_k: departure time from hub k for transportation to destinations
        """
        N = self.set_cities
        return self.model.addVars(N, lb=0, ub=self.par_max_arrival_time_h, vtype = GRB.CONTINUOUS, name="Dc")

    def _add_var_arrival_to_dest(self):
        """A_j: arrival time to destination j
        """
        N = self.set_cities
        return self.model.addVars(N, lb=0, ub=self.par_max_arrival_time_h, vtype = GRB.CONTINUOUS, name="A")


    #########################################
    #                                       #
    #            CONSTRAINTS                #
    #                                       #
    #########################################


    def _add_constraint_assign_node_to_single_hub(self):
        """
        Every node has to be assigned to a single hub

        * Σ_k Z_ik = 1
        * Ɐ i ∈ N
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_assign_orig_hub.sum(i, "*") == 1 for i in N), 
            name = "node to hub assignment")

    def _add_constraint_hub_implication(self):
        """
        Assigning a node to a hub implies the hub was located

        * Z_ik <= Z_kk
        * Ɐ i ∈ N, k ∈ N: i != k
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_assign_orig_hub[i,k] <= self.var_assign_orig_hub[k,k] 
            for i in N for k in N if i != k), 
            name = "node to hub implies hub")

    def _add_constraint_flow_from_source(self):
        """
        The flow from every node goes completely to its hub

        * F_ik = Z_ik * s_i
        * Ɐ i ∈ N, k ∈ N 
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_flow_orig_hub[i,k] == self.var_assign_orig_hub[i,k] * self.par_supply[i] 
             for i in N for k in N),
             name = "flow from source goes to hub")

    def _add_constraint_flow_between_hubs_implies_link(self):
        """
        Flow from origin between hubs implies link between them

        * Y_ikl <= R_kl * s_i
        * Ɐ i ∈ N, Ɐ k ∈ N, l ∈ N: k != l & i != l
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_flow_orig_hub_hub[i,k,l]
             <= 
             self.var_link_hubs[k,l] * self.par_supply[i]
             for i in N for k in N for l in N if (k != l) and (l != i)),
             name = "flow between hubs implies link between them")

    def _add_constraint_conservation_of_flow_at_hub_from_origin(self):
        """
        Conservation flow from commodity (origin) i at hub k: flow in = flow out
        * F_ik + Σ_l (k != l and i != l) Y_ilk  = Σ_l (k != l and i != l) Y_ikl + Σ_j X_ikj
        * Ɐ i ∈ N, Ɐ k ∈ N
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_flow_orig_hub[i,k] +
             gp.quicksum(self.var_flow_orig_hub_hub[i,l,k] for l in N if (i != k) and (l != k))
             == 
             gp.quicksum(self.var_flow_orig_hub_hub[i,k,l] for l in N if (k != l) and (l != i)) + 
             gp.quicksum(self.var_flow_orig_hub_dest[i,k,j] for j in N if i != j)
             for i in N for k in N),
             name = "multicomodity flow conservation at each hub")

    def _add_constraint_flow_node_in_hubs_implies_node_to_hub(self):
        """
        If there is flow from origin i between hubs (k,l) then origin i must have been
        assigned to hub k

        * Y_ikl <= Z_ik * s_i
        * Ɐ i ∈ N, Ɐ k ∈ N, l ∈ N: k != l & l != i
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_flow_orig_hub_hub[i,k,l] <= 
             self.var_assign_orig_hub[i,k] * self.par_supply[i]
             for i in N for k in N for l in N if (k != l) and (l != i)
             ),
             name = "flow from origin in hubs implies origin is assigned to hub"
        )

    def _add_constraint_flow_origin_destination_covered(self):
        """
        The flow from origin i to destination j is covered

        * X_ilj = Z_jl * w_ij
        * Ɐ i ∈ N, Ɐ l ∈ N, j ∈ N: i != j
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_flow_orig_hub_dest[i,l,j] == 
             self.var_assign_orig_hub[j,l] * self.par_flow_box[i,j]
             for i in N for l in N for j in N if i != j
             ),
             name = "flow from origin to destination is covered"
        )

    def _add_constraint_link_between_hubs_implies_hubs(self):
        """
        If there is a link between hubs k and l, then there is a hub in k and in l

        * 2*R_kl <= Z_kk + Z_ll
        * Ɐ k ∈ N, Ɐ l ∈ N: k != l
        """
        N = self.set_cities
        self.model.addConstrs(
            (2*self.var_link_hubs[k,l] <=
             self.var_assign_orig_hub[k,k] + self.var_assign_orig_hub[l,l]
             for k in N for l in N if k != l
             ),
             name = "link between hubs implies hubs"
        )

    def _add_constraint_departure_from_hub_to_hub_waits_customers(self):
        """
        A vehicle departuring from hub k to other hubs must wait cargo from
        its customers (nodes assigned to k) travelling by truck

        * Dh_k >= t_ik * Z_ik
        * Ɐ i ∈ N, Ɐ k ∈ N: i != k
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_departure_from_hub_to_hub[k] >=
             self.par_travel_time_truck_h[i,k] * 
             self.var_assign_orig_hub[i,k]
             for i in N for k in N if i != k
             ),
             name = "departure from hub to hub waits for its customers"
        )
    
    def _add_constraint_departure_from_hub_to_destination_waits_for_hub_trucks(self):
        """
        The departuring time from hub l to its customers (DHC_l) must
        wait for incoming cargo from other hubs k Dh_k to consolidate cargo.
        
        * Original non-linear: Dc_l >= (Dh_k + t_kl*α) * R_kl
        * Linearized: Dc_l >= Dh_k + t_kl*α - M(1 - R_kl)
        * Where M is a sufficient big number, in this case: M = β (max delivery time)
        * Ɐ k ∈ N, Ɐ l ∈ N
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_departure_from_hub_to_dest[l] >=
             self.var_departure_from_hub_to_hub[k] + 
             self.par_travel_time_truck_h[k,l] * self.par_economy_scale_factor -
             self.par_max_arrival_time_h*(1 - self.var_link_hubs[k,l])
             for k in N for l in N
             ),
             name = "departure from hub to destination waits for hub truck"
        )
    
    def _add_constraint_departure_from_hub_to_destination_waits_for_customer_trucks(self):
        """
        The departuring time from hub l to its customers j (Dc_l) must
        wait for incoming truck cargo from its customers to consolidate. This is a 
        special case when there is only a single hub in the network.
        
        * Dc_l >= t_jl * Z_jl
        * Ɐ j ∈ N, Ɐ l ∈ N
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_departure_from_hub_to_dest[l] >=
             self.par_travel_time_truck_h[j,l] *
             self.var_assign_orig_hub[j,l]
             for j in N for l in N
             ),
             name = "departure from hub to destination waits for truck customers"
        )

    def _add_constraint_latest_arrival_time_to_destination(self):
        """
        The latest arrival time to destination j depends on the latest truck
        arriving from hub l

        * Original non-linear: A_j >= (Dc_l + t_lj) * Z_jl
        * Linearized: A_j >= Dc_l + t_lj - M(1 - Z_jl)
        * Where M is a sufficient big number, in this case: M = β (max delivery time)
        * Ɐ j ∈ N, Ɐ l ∈ N
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_arrival_to_dest[j] >=
             self.var_departure_from_hub_to_dest[l] + 
             self.par_travel_time_truck_h[l,j] - 
             self.par_max_arrival_time_h*(1 - self.var_assign_orig_hub[j,l])
             for j in N for l in N
             ),
             name = "departure from hub to destination waits for hub to hub trucks"
        )

    def _add_constraint_latest_arrival_to_dest_within_limit(self):
        """
        The latest arrival time to destination j must be within the standard

        * A_j <= β
        * Ɐ j ∈ N
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_arrival_to_dest[j] <= self.par_max_arrival_time_h
             for j in N
             ),
             name = "max arrival time within standard"
        )

    def _add_contraint_heuristic_non_hubs(self):
        """We limit here the possibility of some nodes to be hubs by the precomputed parameter
        heuristic_non_hubs_nodes
        - Z_ii = 0 
        - Ɐ i ∈ Discarded Hubs
        """
        self.model.addConstrs(
            (self.var_assign_orig_hub[i,i] == 0
             for i in self.heuristic_non_hubs_nodes
             ),
             name = "heuristic discarded hubs"
        )


    #########################################
    #                                       #
    #                  CUTS                 #
    #                                       #
    #########################################
    
    # def _add_constraint_eliminate_slowest_combinations(self):
    #     """
    #     There are some combinations of routes that are infeasible because
    #     they are too slow for the given time standard β. These are such
    #     routes i -> k -> l -> j, where k and l are hubs. The condition is:

    #     * if tt_ik + tp_kl + tt_lj > β
    #     * then: Z_ik + R_kl + Z_jl <= 2
    #     * Ɐ i ∈ N, Ɐ k ∈ N, Ɐ l ∈ N, Ɐ j ∈ N: k != l, i != j
    #     """
    #     N = self.set_cities
    #     tt = self.par_travel_time_truck_h
    #     tp = self.par_travel_time_plane_h
    #     beta = seleconomy_scale_factor
    #     infeasible_combinations = []
    #     print(f"\nApplying cut: detecting tt_ik + tp_kl + tt_lj > β")
    #     for i in N:
    #         for k in N:
    #             for l in N:
    #                 for j in N:
    #                     if (k != l) and (i != j):
    #                         time_iklj = tt[i,k] + tp[k,l] + tp[l,j]
    #                         if time_iklj > beta:
    #                             infeasible_combinations.append((i,k,l,j))
    #                             print(f"\ttime through {i,k,l,j} = {time_iklj} > {beta}")
    #     print(f"\tTotal number of this cut: {len(infeasible_combinations)}")

    #     self.model.addConstrs(
    #         (self.var_assign_orig_hub[i,k] +
    #          self.var_link_hubs[k,l] + 
    #          self.var_assign_orig_hub[j,l]
    #          <=
    #          2
    #          for (i,k,l,j) in infeasible_combinations
    #          ),
    #          name = "cut routes taking longer than β"
    #     )

    



    #########################################
    #                                       #
    #          OBJECTIVE FUNCTION           #
    #                                       #
    #########################################

    def _add_var_OF_location_cost(self):
        return self.model.addVar(lb=0, vtype = GRB.CONTINUOUS, name="OF_location_cost")
    
    def _add_constraint_OF_location_cost(self):
        """The location cost of the hubs

        * Σ_i Z_ii * h_i
        """
        N = self.set_cities
        self.model.addConstr(
            self.var_OF_location_cost == 
            gp.quicksum(
                self.var_assign_orig_hub[i,i] * self.par_fixed_hub_cost[i]
                for i in N
            ),
            name = "OF_location_cost"
        )

    def _add_var_OF_transport_collection_cost(self):
        return self.model.addVar(lb=0, vtype = GRB.CONTINUOUS, name="OF_transport_collection_cost")

    def _add_constraint_OF_transport_collection_cost(self):
        """The transportation cost of collection activities

        * Σ_i Σ_k F_ik * d_ik * c
        """
        N = self.set_cities
        self.model.addConstr(
            self.var_OF_transport_collection_cost == 
            gp.quicksum(
                self.var_flow_orig_hub[i,k] * self.par_distance_km[i,k] * self.par_unit_cost_per_km_flow
                for i in N for k in N
            ),
            name = "OF_transport_collection_cost"
        )

    def _add_var_OF_transport_transfer_cost(self):
        return self.model.addVar(lb=0, vtype = GRB.CONTINUOUS, name="OF_transport_transfer_cost")

    def _add_constraint_OF_transport_transfer_cost(self):
        """The transportation cost of transfer activities

        * Σ_i Σ_k Σ_l (k != l & l != i) Y_ikl * d_kl * c * α
        """
        N = self.set_cities
        self.model.addConstr(
            self.var_OF_transport_transfer_cost == 
            gp.quicksum(
                self.var_flow_orig_hub_hub[i,k,l] * 
                self.par_distance_km[k,l] * 
                self.par_unit_cost_per_km_flow * 
                self.par_economy_scale_factor
                for i in N for k in N for l in N if (k != l) and (l != i)
            ),
            name = "OF_transport_transfer_cost"
        )

    def _add_var_OF_transport_distribution_cost(self):
        return self.model.addVar(lb=0, vtype = GRB.CONTINUOUS, name="OF_transport_distribution_cost")
    
    def _add_constraint_OF_transport_distribution_cost(self):
        """The transportation cost of distribution activities

        * Σ_i Σ_l Σ_j X_ilj * d_lj * c
        """
        N = self.set_cities
        self.model.addConstr(
            self.var_OF_transport_distribution_cost == 
            gp.quicksum(
                self.var_flow_orig_hub_dest[i,l,j] * self.par_distance_km[l,j] * self.par_unit_cost_per_km_flow
                for i in N for l in N for j in N if i != j
            ),
            name = "OF_transport_distribution_cost"
        )

    def _add_var_OF_links_cost(self):
        return self.model.addVar(lb=0, vtype = GRB.CONTINUOUS, name="OF_links_cost")

    def _add_constraint_OF_links_cost(self):
        """The fixed cost of stablishing links between nodes

        * Σ_i Σ_k Z_ik * f_ik + Σ_k Σ_l R_kl * f_kl
        """
        N = self.set_cities
        self.model.addConstr(
            self.var_OF_links_cost == 
            gp.quicksum(
                self.var_assign_orig_hub[i,k] * self.par_fixed_link_cost[i,k]
                for i in N for k in N
            ) + 
            gp.quicksum(
                self.var_link_hubs[k,l] * self.par_fixed_link_cost[k,l]
                for k in N for l in N
            ),
            name = "OF_links_cost"
        )
    
    def _set_objective_function(self):
        self.model.setObjective(
            self.var_OF_links_cost +
            self.var_OF_location_cost + 
            self.var_OF_transport_collection_cost + 
            self.var_OF_transport_distribution_cost + 
            self.var_OF_transport_transfer_cost,
            sense = GRB.MINIMIZE
        )

    def _set_start_values(self):
        """Here we set some initial values to speed up finding a feasible solution
        """
        # For each potential node that can be a hub, set it as a hub
        # For each non-hub node, assign it to the closest hub
        
        if self.heuristic_non_hubs_nodes is not None:
            print(f"\nApplying heuristic start values")
            hub_ids = set(self.set_cities.copy()) - self.heuristic_non_hubs_nodes
            for k in hub_ids:
                self.var_assign_orig_hub[k,k].Start = 1.0
            print(f"Hubs are: {list(hub_ids)}")
            
            # Assign these cities to the closest hub
            print("Assignment to hubs:")
            for i in self.heuristic_non_hubs_nodes:
                time_hub = [(self.par_travel_time_truck_h[(i,k)], k) for k in hub_ids]
                time_hub.sort()
                k = time_hub[0][1]
                self.var_assign_orig_hub[i,k].Start = 1.0
                print(f"\tNode {i} assigned to Hub {k}")


    def solve(self):
        print("\n\n\nSOLVING...\n")
        self.model.optimize()

        if self.model.status == GRB.INFEASIBLE:
            print(f"\n\nModel infeasable, computing IIS:")
            self.model.computeIIS()
            report = self.instance_name + "_IIS.ilp"
            self.model.write(report)
            print(f"IIS report written in {report}")

    def save_solution(self):
        """Saves the variables with values != 0 in a json file
        """
        tolerance = float("1.0e-10")
        def remove_0_values(d: dict()) -> dict():
            return { str(k): v for (k,v) in d.items() if v > tolerance}

        solution = dict()

        solution["MIPGap"] = self.model.getAttr("MIPGap")
        
        # Objective function value
        solution["OF_val"] = self.model.ObjVal

        # Time needed to solve problem
        solution["time_s"] = self.model.Runtime

        # Append objective function components
        solution['OF'] = {}
        solution['OF']['location_cost'] = self.var_OF_location_cost.X
        solution['OF']['transport_collection_cost'] = self.var_OF_transport_collection_cost.X
        solution['OF']['transport_transfer_cost'] = self.var_OF_transport_transfer_cost.X
        solution['OF']['transport_distribution_cost'] = self.var_OF_transport_distribution_cost.X
        solution['OF']['links_cost'] = self.var_OF_links_cost.X

        # Append variables
        solution['variables'] = {}
        solution['variables']['F'] = remove_0_values(self.model.getAttr('X', self.var_flow_orig_hub))
        solution['variables']['Z'] = remove_0_values(self.model.getAttr('X', self.var_assign_orig_hub))
        solution['variables']['Y'] = remove_0_values(self.model.getAttr('X', self.var_flow_orig_hub_hub))
        solution['variables']['X'] = remove_0_values(self.model.getAttr('X', self.var_flow_orig_hub_dest))
        solution['variables']['R'] = remove_0_values(self.model.getAttr('X', self.var_link_hubs))
        solution['variables']['Dh'] = remove_0_values(self.model.getAttr('X', self.var_departure_from_hub_to_hub))
        solution['variables']['Dc'] = remove_0_values(self.model.getAttr('X', self.var_departure_from_hub_to_dest))
        solution['variables']['A'] = remove_0_values(self.model.getAttr('X', self.var_arrival_to_dest))

        utils.save_solution(solution, file_name=self.instance_name)

    def load_solution(self):
        """Loads the values of the variables in a saved solution
        """        
        return utils.load_solution(self.instance_name)

    def plot_solution(self):
        solution = self.load_solution()
        
        # Input for plotting
        hubs_ids = set()
        origin_hub_flow_collection = solution["variables"]["F"]
        origin_hub_hub_flow_transfer = solution["variables"]["Y"]
        origin_hub_destination_flow_distribution = solution["variables"]["X"]

        # Getting hubs_ids
        Z = solution["variables"]["Z"]
        for (i,j), _ in Z.items():
            if i == j:
                hubs_ids.add(j)
        
        # Now we can plot
        potential_hubs = self.set_cities - self.heuristic_non_hubs_nodes
        plotter.plot_map(cities_data = self.cities_data,
                         cities_considered = self.set_cities,
                         hubs_ids=hubs_ids,
                         collection=origin_hub_flow_collection,
                         transfer=origin_hub_hub_flow_transfer,
                         distribution=origin_hub_destination_flow_distribution,
                         beta=self.par_max_arrival_time_h,
                         alpha=self.par_economy_scale_factor,
                         potential_hubs=potential_hubs,
                         top_k_cities_for_hub=self.top_k_cities_for_hub,
                         size_proportional_to_flow=True,
                         draw_city_names=False
                         )



if __name__ == "__main__":
    from dataloader import load_data
    cities_data = load_data()
    
    problem = USAPSCT(cities_data,
                      max_nodes=81,
                      max_arrival_time_h=29,
                      economy_of_scale_factor=0.8,
                      top_k_cities_for_hub=20
                      )
    problem.solve()
    problem.save_solution()
    problem.plot_solution()

    # Experiment 1: change beta
    # betas = [20, 24, 28, 99] 
    # for b in betas:
    #     problem = USAPSCT(cities_data,
    #                     max_nodes=81,
    #                     number_hubs=None,
    #                    top_k_nodes_hubs
    #                     hub_to_hub_cost_factor=8, 
    #                     min_city_supply_ranking_to_be_hub=41,
    #                     min_latitude_to_be_hub=37,
    #                     max_latitude_to_be_hub=41.02,
    #                     min_longitude_to_be_hub=28,
    #                     max_longitude_to_be_hub=41)
    #     problem.solve()
    #     problem.save_solution()
    #     problem.plot_solution()
    
    # Experiment 2: change alpha
    # alphas = [2, 4, 8, 16] 
    # for a in alphas:
    #     problem = USAPSCT(cities_data,
    #                     max_nodes=81,
    #                     number_hubs=None,
    #                     max_time_h=24,
    #                     hub_to_hub_cost_factor=a, 
    #                     min_city_supply_ranking_to_be_hub=41,
    #                     min_latitude_to_be_hub=37,
    #                     max_latitude_to_be_hub=41.02,
    #                     min_longitude_to_be_hub=28,
    #                     max_longitude_to_be_hub=41)
    #     problem.solve()
    #     problem.save_solution()
    #     problem.plot_solution()
    
    # Experiment 3: look at heuristic effect (need last 2)
    # params = [
    #     (60, 38.0, 40.0, 32.0, 40.0),
    #     (60, 37.5, 40.5, 31.5, 40.5),
    #     (60, 37.0, 41.0, 31.0, 41.0),
    #     (60, 36.5, 41.5, 30.5, 41.5),
    #     (60, 36.0, 42.0, 30.0, 42.0),
    #     (60, 35.5, 42.0, 29.5, 42.5)
    # ] 
    # for p in params:
    #     min_ranking, min_lat, max_lat, min_long, max_long = p
    #     problem = USAPSCT(cities_data,
    #                     max_nodes=81,
    #                     number_hubs=None,
    #                     max_time_h=24,
    #                     hub_to_hub_cost_factor=8, 
    #                     min_city_supply_ranking_to_be_hub=min_ranking,
    #                     economy_of_scale_factor=min_lat,
    #                     max_latitude_to_be_hub=max_lat,
    #                     min_longitude_to_be_hub=min_long,
    #                     max_longitude_to_be_hub=max_long)
    #     problem.solve()
    #     problem.save_solution()
    #     problem.plot_solution()
