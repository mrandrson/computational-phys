#include <stdio.h>
#include <stdlib.h>

typedef struct {
    float x;
    float y;
} Point;

typedef struct QuadTree {
    Point* data;               
    int data_count;            
    float min_x, min_y;        
    float max_x, max_y;
    int max_iterations;        
    struct QuadTree* children[4]; 
} QuadTree;

QuadTree* createQuadTree(Point* data, int data_count, float min_x, float min_y, float max_x, float max_y, int max_iterations);
void subdivideQuadTree(QuadTree* qt);
void freeQuadTree(QuadTree* qt);

QuadTree* createQuadTree(Point* data, int data_count, float min_x, float min_y, float max_x, float max_y, int max_iterations) {
    QuadTree* qt = (QuadTree*)malloc(sizeof(QuadTree));
    qt->data = (Point*)malloc(data_count * sizeof(Point));
    for (int i = 0; i < data_count; i++) {
        qt->data[i] = data[i];
    }
    qt->data_count = data_count;
    qt->min_x = min_x;
    qt->min_y = min_y;
    qt->max_x = max_x;
    qt->max_y = max_y;
    qt->max_iterations = max_iterations;

    for (int i = 0; i < 4; i++) {
        qt->children[i] = NULL;
    }

    if (max_iterations > 0 && data_count > 0) {
        subdivideQuadTree(qt);
    }

    return qt;
}

/*
void subdivideQuadTree(QuadTree* qt) {
    if (qt->max_iterations <= 0 || qt->data_count == 0) return;

    float mid_x = (qt->min_x + qt->max_x) / 2;
    float mid_y = (qt->min_y + qt->max_y) / 2;

    Point* q1_data = malloc(qt->data_count * sizeof(Point));
    Point* q2_data = malloc(qt->data_count * sizeof(Point));
    Point* q3_data = malloc(qt->data_count * sizeof(Point));
    Point* q4_data = malloc(qt->data_count * sizeof(Point));

    int q1_count = 0, q2_count = 0, q3_count = 0, q4_count = 0;

    for (int i = 0; i < qt->data_count; i++) {
        Point p = qt->data[i];
        if (p.x < mid_x && p.y < mid_y) q1_data[q1_count++] = p;
        else if (p.x < mid_x && p.y >= mid_y) q2_data[q2_count++] = p;
        else if (p.x >= mid_x && p.y < mid_y) q3_data[q3_count++] = p;
        else if (p.x >= mid_x && p.y >= mid_y) q4_data[q4_count++] = p;
    }

    qt->children[0] = createQuadTree(q1_data, q1_count, qt->min_x, qt->min_y, mid_x, mid_y, qt->max_iterations - 1);
    qt->children[1] = createQuadTree(q2_data, q2_count, qt->min_x, mid_y, mid_x, qt->max_y, qt->max_iterations - 1);
    qt->children[2] = createQuadTree(q3_data, q3_count, mid_x, qt->min_y, qt->max_x, mid_y, qt->max_iterations - 1);
    qt->children[3] = createQuadTree(q4_data, q4_count, mid_x, mid_y, qt->max_x, qt->max_y, qt->max_iterations - 1);

    free(q1_data);
    free(q2_data);
    free(q3_data);
    free(q4_data);
}
*/

void subdivideQuadTree(QuadTree* qt) {
    if (qt->max_iterations <= 0 || qt->data_count == 0) return;

    float mid_x = (qt->min_x + qt->max_x) / 2;
    float mid_y = (qt->min_y + qt->max_y) / 2;

    // Allocate memory for the four quadrants
    Point* q1_data = malloc(qt->data_count * sizeof(Point));
    Point* q2_data = malloc(qt->data_count * sizeof(Point));
    Point* q3_data = malloc(qt->data_count * sizeof(Point));
    Point* q4_data = malloc(qt->data_count * sizeof(Point));

    int q1_count = 0, q2_count = 0, q3_count = 0, q4_count = 0;

    // Distribute points into quadrants
    for (int i = 0; i < qt->data_count; i++) {
        Point p = qt->data[i];
        if (p.x < mid_x && p.y < mid_y) q1_data[q1_count++] = p;
        else if (p.x < mid_x && p.y >= mid_y) q2_data[q2_count++] = p;
        else if (p.x >= mid_x && p.y < mid_y) q3_data[q3_count++] = p;
        else if (p.x >= mid_x && p.y >= mid_y) q4_data[q4_count++] = p;
    }

    // Debug prints
    /*
    printf("Subdividing node with bounds [%f, %f] to [%f, %f]\n", qt->min_x, qt->min_y, qt->max_x, qt->max_y);
    printf("Quadrant 1: %d points\n", q1_count);
    printf("Quadrant 2: %d points\n", q2_count);
    printf("Quadrant 3: %d points\n", q3_count);
    printf("Quadrant 4: %d points\n", q4_count);
    */

    // Create child nodes with copied data
    qt->children[0] = createQuadTree(q1_data, q1_count, qt->min_x, qt->min_y, mid_x, mid_y, qt->max_iterations - 1);
    qt->children[1] = createQuadTree(q2_data, q2_count, qt->min_x, mid_y, mid_x, qt->max_y, qt->max_iterations - 1);
    qt->children[2] = createQuadTree(q3_data, q3_count, mid_x, qt->min_y, qt->max_x, mid_y, qt->max_iterations - 1);
    qt->children[3] = createQuadTree(q4_data, q4_count, mid_x, mid_y, qt->max_x, qt->max_y, qt->max_iterations - 1);

    // Free temporary arrays (ownership transferred to children)
    free(q1_data);
    free(q2_data);
    free(q3_data);
    free(q4_data);
}


void freeQuadTree(QuadTree* qt) {
    if (qt == NULL) return;

    for (int i = 0; i < 4; i++) {
        freeQuadTree(qt->children[i]);
    }

    free(qt->data);
    free(qt);
}

