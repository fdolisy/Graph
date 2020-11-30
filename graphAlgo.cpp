/**
* CS 3345 HON
* Project 3, GRAPHS
* Used compiler g++ 4.8.5 - 11 with compile command : g++ -std=c++11 -o CMLQP3 CMLQP3.cpp
* date 11 / 20 / 20
*/
//NOTE: on average takes about 5 seconds for large test files 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits.h>
#include <math.h>
#include <algorithm>
//#include <ios>
//#include <time.h>

//a function used to round numbers to 2 decimal places
float roundx(float val)
{
    float value = (int)(val * 100 + .5);
    return (float) value / 100;
}
//class for edges
class Edge
{
    public:
        int x;
        int y;
        double w;
    //creating my edge
    Edge(int x, int y, double w)
    {
        this->x = x;
        this->y = y;
        this->w = w;
    }
    //default constructor
    Edge()
    {
        x = 0;
        y = 0;
        w = 0;
    }

};
//used to compare edges in Kruskals(std::sort)
bool operator<(Edge const& a, Edge const& b)
{
    if (a.w < b.w)
    {
        return true;
    }
    else if (a.w > b.w)
    {
        return false;
    }
    else
    {
        return a.x < b.x;
    }
}

//class for disjoint sets(used for Kruskals)
class DisjointSet
{
    public:
        int* rank;
        int* parent;
        int n;

        DisjointSet(int n);
        int find(int u);
        void graphUnion(int x, int y);
   
};
//constructor
DisjointSet::DisjointSet(int n)
{

    this->n = n;
    parent = new int[n + 1];
    rank = new int[n + 1];


    for (int i = 0; i <= n; i++)
    {
        rank[i] = 0;
        parent[i] = i;
    }
}
//find the parent of the input
int DisjointSet::find(int u)
{
    if (u != parent[u])
        parent[u] = find(parent[u]);
    return parent[u];
}
// Union by rank 
void DisjointSet::graphUnion(int x, int y)
{
    x = find(x);
    y = find(y);
    if (rank[x] > rank[y])
    {
        parent[y] = x;
    }
    else
    {
        parent[x] = y;
    }
       

    if (rank[x] == rank[y])
    {
        rank[y]++;
    }
}

//the graph class
class Graph
{
    public:

        //variables used throughout
        double** adjMatrix;
        int numVertices;
        int** dist;
        int** path;
        int edgeCount;
        Edge* e;
    
        //basic utility functions
        Graph(int);
        ~Graph();
        void setEdgeCount(int);
        void addEdge(int, int, double);
        void display();
        void transform();

        //key calculating algorithms
        float kruskal();
        int floyd();
        void printPath(int u, int v, int);
        int color();
       
};

//constructor of graph
Graph::Graph(int numVertices)
{
    this->numVertices = numVertices;
    this->edgeCount = 0;

    e = new Edge[(numVertices * numVertices - 1) / 2];
    dist = new int* [numVertices];
    path = new int* [numVertices];

    adjMatrix = new double* [numVertices];
    for (int i = 0; i < numVertices; i++)
    {
        dist[i] = new int[numVertices];
        adjMatrix[i] = new double[numVertices];
        path[i] = new int[numVertices];
        for (int j = 0; j < numVertices; j++)
        {
            adjMatrix[i][j] = 0;
            dist[i][j] = 0;
            path[i][j] = 0;
        }
    }
   
}
//destructor of graph
Graph::~Graph()
{
    for (int i = 0; i < numVertices; i++)
    {
        delete[] adjMatrix[i];
        delete[] dist[i];
        delete[] path[i];
    }
    
    delete[] dist;
    delete[] path;
    delete[] adjMatrix;
    delete[] e;
}
//sets the edge count (used for the edge array)
void Graph::setEdgeCount(int x)
{
    edgeCount = x;

}
//adds an edge to the graphwith its weight
void Graph::addEdge(int x, int y, double val)
{
   
    adjMatrix[x][y] = adjMatrix[y][x] = val;
    //adjMatrix[y][x] = val;
    Edge* z = new Edge(x, y, val);
    e[edgeCount] = *z;

    edgeCount++;

}
//display the adjcency matrix of the graph (for debugging purposes)
void Graph::display()
{
    for (int i = 0; i < numVertices; i++)
    {
        for (int j = 0; j < numVertices; j++)
        {
            std::cout << adjMatrix[i][j] << " ";
        }

        // Newline for new row 
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for (int i = 0; i < numVertices; i++)
    {
        for (int j = 0; j < numVertices; j++)
        {
            std::cout << dist[i][j] << " ";
        }

        // Newline for new row 
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for (int i = 0; i < numVertices; i++)
    {
        for (int j = 0; j < numVertices; j++)
        {
            std::cout << path[i][j] << " ";
        }

        // Newline for new row 
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
//transform the weights of the graphs to 1's and zero's
void Graph::transform()
{
    for (int i = 0; i < numVertices; i++)
    {
        for (int j = 0; j < numVertices; j++)
        {
            if (adjMatrix[i][j] != 0)
            {
                if (adjMatrix[i][j] == INT_MAX)
                {
                    adjMatrix[i][j] = 0;
                    dist[i][j] = adjMatrix[i][j];
                    path[i][j] = j;
                }
                else
                {
                    adjMatrix[i][j] = 1;
                    dist[i][j] = adjMatrix[i][j];
                    path[i][j] = j;
                }
            }
            if (i == j)
            {
                path[i][j] = i;
            }

        }
    }
}
// return the total weight of the MST and print the individual weights
float Graph::kruskal()
{
    double minCost = 0, min = 0; 
    int edgeC = 0;
    
    //sort the edges according to weight
    std::sort(e, e + edgeCount);
      
    // create disjoint sets 
    DisjointSet ds(numVertices);

    while (edgeC < numVertices - 1)
    {

        for (int i = 0; i < edgeCount; i++)
        {
            int a = e[i].x;
            int b = e[i].y;

            int aSet = ds.find(a);
            int bSet = ds.find(b);

            //no cycle will be created
            if (aSet != bSet)
            {

                min = e[i].w;

                std::cout << a + 1 << " " << b + 1 << " " << std::fixed << std::setprecision(2) << roundx(min) << std::endl;

                minCost += min;
                ds.graphUnion(aSet, bSet);
                edgeC++;
            }
        }
    }
    return roundx(minCost);
}
//calculates the diameter and finds calculate all paths
int Graph::floyd()
{
    int n = numVertices;
    int di = 0;

    //basic implementation of kruskals for unweighted graph
    for (int k = 0; k < n; ++k)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                if ((dist[i][k] != 0) && (dist[k][j] != 0) && (i != j))
                {
                    if ((dist[i][k] + dist[k][j] < dist[i][j]) || (dist[i][j] == 0))
                    {
                        dist[i][j] = dist[j][i] = dist[i][k] + dist[k][j];
                        path[i][j] = path[j][i] = path[i][k];
                    }
                }
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            //graph is symettric so we only need to visit upper triangular to find max
            if (i < j && i != 0)
            {
                break;
            }
            else
            {
                if (dist[i][j] > di)
                {
                    di = dist[i][j];
                }
                if (i == 0 && i != j)
                {
                    printPath(i, j, dist[i][j]);
                }
            }
        }
    }
    return di;
}
//prints out the shortest path
void Graph::printPath(int u, int v, int d)
{
    std::cout << u + 1 << " ";
    while (u != v)
    {
        u = path[u][v];
        std::cout << u + 1 << " ";
    }
    std::cout << d << std::endl;

}
//colors the graph and returns the number of colors used
int Graph::color()
{
    int numColors = 0;
    int* color = new int[numVertices];
    bool* used = new bool[numVertices];
    color[0] = 0;
    
    //initialize all other vertices as empty
    for (int i = 1; i < numVertices; i++)
    {
        color[i] = -1;
    }

    //initialize all colors as unused
    for (int i = 0; i < numVertices; i++)
    {
        used[i] = false;
    }

    //for all other numVertices - 1 vertices
    for (int i = 1; i < numVertices; i++)
    {
        for (int j = 0; j < numVertices; j++) 
        {
            if (adjMatrix[i][j]) {
                if (color[j] != -1)
                    used[color[j]] = true;
            }
        }

        int a;
        for (a = 0; a < numVertices; a++)
        {
            if (!used[a])
            {
                break;
            }
        }

        color[i] = a;    //assign found color in the list

        for (int k = 0; k < numVertices; k++) 
        {
            if (adjMatrix[i][k])
            {
                if (color[k] != -1)
                {
                    used[color[k]] = false;
                }
            }
        }
    }

    //get the number of colors used
    for (int u = 0; u < numVertices; u++)
    {
        int j = 0;
        for (j = 0; j < u; j++)
            if (color[u] == color[j])
            {
                break;
            }

   
        if (u == j)
        {
            numColors++;
        }
    }

    return numColors;
}

int main()
{
    //to measure time running
    /*
    time_t start, end;
    time(&start);
    std::ios_base::sync_with_stdio(false);
    */

    std::ifstream myfile;
    int num;
    int edgeCout = 0;
    double x = 0, y = 0, radius = 0, sqrtr = 0;

    myfile.open("GraphData.txt");
    myfile >> num;
    double* xar = new double[num];
    double* yar = new double[num];
    Graph g(num);

    //save all x and y values to be compared later
    for (int i = 0; i < num; i++)
    {
        myfile >> x >> y;
        xar[i] = x;
        yar[i] = y;
    }

    //get the radius
    myfile >> radius;

    //compare all vertices to one another
    for (int i = 0; i < num; i++) 
    {
        for (int k = 1 + i; k < num; k++) 
        {
            sqrtr = sqrt((pow(xar[i] - xar[k], 2) + (pow(yar[i] - yar[k], 2))));
            //if the distance is less than the radius, create an edge
            if (sqrtr <= radius)
            {
                g.addEdge(i, k, sqrtr);
                edgeCout++;

            }
            else
            {
                g.addEdge(i, k, INT_MAX);
                edgeCout++;
            }
        }
    }

    g.setEdgeCount(edgeCout);

    //what gets outputted
    std::cout << g.kruskal() << std::endl;
    g.transform();
    std::cout << g.floyd() << std::endl;
    std::cout << g.color() << std::endl;

    myfile.close();

    /*
    time(&end);
    double time_taken = double(end - start);
    std::cout << "Time taken by program is : " << std::fixed
        << time_taken << std::setprecision(5);
    std::cout << " sec " << std::endl;
    */
}
