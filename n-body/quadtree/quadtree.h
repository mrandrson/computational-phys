#ifndef QUADTREE_H
#define QUADTREE_H

// Struct for a 2D point
typedef struct {
    float x;
    float y;
} Point;

// Struct for a QuadTree node
typedef struct QuadTree {
    Point* data;               // Array of points in this node
    int data_count;            // Number of points in this node
    float min_x, min_y;        // Bounds
    float max_x, max_y;
    int max_iterations;        // Remaining iterations
    struct QuadTree* children[4]; // Pointers to child nodes (Q1, Q2, Q3, Q4)
} QuadTree;

// Function to create a QuadTree node
QuadTree* createQuadTree(Point* data, int data_count, float min_x, float min_y, float max_x, float max_y, int max_iterations);

// Function to subdivide a QuadTree node
void subdivideQuadTree(QuadTree* qt);

// Function to free the QuadTree
void freeQuadTree(QuadTree* qt);

#endif // QUADTREE_H

