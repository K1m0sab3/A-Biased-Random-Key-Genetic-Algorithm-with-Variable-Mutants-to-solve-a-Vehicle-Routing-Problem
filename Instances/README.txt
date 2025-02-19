Instance sets used in Festa P. et al., "A Biased Random-Key Genetic Algorithm with Variable Mutants to solve a Vehicle Routing Problem".

We used the same problem instances analyzed by Macrina et al. in their work "A Variable Neighborhood Search for the
Vehicle Routing Problem with Occasional Drivers and Time Windows". However, we modified the data format to better suit
our needs. To facilitate the computational analysis and experimentation, we grouped the instances into four distinct
classes, ordered by increasing size. Each class contains three different network typologies: clustered-type, random-type,
and mixed-type. The first class (Test) consists of a total of 36 instances, with 12 instances each for 5, 10,and 15 customers,
respectively. Each of the other three classes (Small, Medium, and Large) consists of 15 instances with 25, 50,
and 100 customers, respectively. The instance characteristics are shown in the paper.

Each instance is a .dat file formatted as follows:

1. num nodes, num drivers, num Occasional Drivers (OD), capacities of occasional drivers, capacities of classical drivers (C)

2. An upper triangular portion of matrix (TRIU) of distances as following (for n customers and m OD)

        dist(C1, C2) dist(C1, C3) ... dist(C1, Cn) ||| dist(C1, OD1) ... dist(C1, ODm) dist(C1, Depot)
        dist(C2, C3) ... dist(C2, Cn) ||| dist(C2, OD1) ... dist(C2, ODm) dist(C2, Depot)
        ...
        dist(Cn-1, Cn) ||| dist(Cn-1, OD1) ... dist(Cn-1, ODm) dist(Cn-1, Depot)
        ||| dist(Cn, OD1) ... dist(Cn, ODm) dist(Cn, Depot)

Note: The entire matrix is not actually a TRIU. To identify the TRIU structure, you should focus only on the left section of the matrix, marked by the three vertical bars |||.

This TRIU section represents the distances between customers, while the rest of the matrix contains distances:
- Between customers and ODs
- Between customers and the depot.

3. A (n+m)X2 matrix of time windows for the n customers and the m occasional drivers

        lC1 uC1 (lower-time windows-customer-1 upper-time windows-customer-1)
        ... ...     ...
        lCn  uCn (lower-time windows-customer-n upper-time windows-customer-n)
        lOD1 uOD1 (lower-time windows-occasional-driver-1 upper-time windows-occasional-driver-1)
        ...  ...    ...
        lODm uODm (lower-time windows-occasional-driver-m upper-time windows-occasional-driver-m)

4. Distances from the depot

5. Customer requests
