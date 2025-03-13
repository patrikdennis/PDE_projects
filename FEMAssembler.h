#ifndef FEM_ASSEMBLER_H
#define FEM_ASSEMBLER_H

#include <vector>
#include <array>
#include <functional>

// Structure to represent a 2D point.
struct Node {
    double x, y;
};

class FEMAssembler {
private:
    // Mesh data: nodes, elements (triangles), and boundary edges.
    std::vector<Node> nodes;                         // Each node: x and y coordinates.
    std::vector<std::array<int, 3>> elements;          // Each element (triangle): indices of 3 nodes.
    std::vector<std::array<int, 2>> boundaries;        // Each boundary edge: indices of 2 nodes.

public:
    // Constructor (assumes indices are 0-based).
    FEMAssembler(const std::vector<Node>& nodes,
                 const std::vector<std::array<int, 3>>& elements,
                 const std::vector<std::array<int, 2>>& boundaries);

    // Assemble the boundary stiffness matrix.
    std::vector<std::vector<double>> assembleBoundaryMatrix(double c);

    // Assemble the boundary load vector.
    std::vector<double> assembleBoundaryVector(double c, double u0);

    // Assemble the stiffness matrix for constant thermal conductivity.
    std::vector<std::vector<double>> assembleStiffnessMatrix(double a);

    // Assemble the stiffness matrix using nodal variable conductivity.
    std::vector<std::vector<double>> assembleStiffnessMatrix2(const std::vector<double>& aCoeff);

    // Assemble the stiffness matrix using quadrature (centroid evaluation).
    std::vector<std::vector<double>> assembleStiffnessMatrixQuadrature(const std::vector<double>& aCoeff);

    // Assemble the load vector using a quadrature rule.
    std::vector<double> assembleLoadVectorQuadrature(std::function<double(double, double)> f_handle);
};

#endif // FEM_ASSEMBLER_H

