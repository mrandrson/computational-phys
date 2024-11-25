#include <stdio.h>
#include "quadtree.h"

int main() {
    Point data[] = {
        {1, 1}, {2, 3}, {4, 2}, {3, 5}, {6, 7}, {7, 8}, {14, 15}
    };
    int data_count = sizeof(data) / sizeof(data[0]);

    QuadTree* qt = createQuadTree(data, data_count, 0, 0, 10, 10, 3);

    printf("Root Node:\n");
    printf("Bounds: [%f, %f] to [%f, %f]\n", qt->min_x, qt->min_y, qt->max_x, qt->max_y);
    printf("Data Count: %d\n", qt->data_count);

    freeQuadTree(qt);

    return 0;
}

