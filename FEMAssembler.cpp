#include <vector>
#include <array>
#include <cmath>
#include <functional>

// hold 2D point coordinates
struct Node {
    double x, y;
};

class FEMAssembler {
private:
    // Mesh data: nodes, elements (triangles), and boundary edges.
    std::vector<Node> nodes;                         // Each node has x and y coordinates.
    std::vector<std::array<int, 3>> elements;          // Each element (triangle): indices of 3 nodes.
    std::vector<std::array<int, 2>> boundaries;        // Each boundary edge: indices of 2 nodes.

public:
    // Constructor: assume indices are 0-based.
    FEMAssembler(const std::vector<Node>& nodes,
                 const std::vector<std::array<int, 3>>& elements,
                 const std::vector<std::array<int, 2>>& boundaries)
        : nodes(nodes), elements(elements), boundaries(boundaries) {}

    // Assemble the boundary stiffness matrix B.
    // c is a scaling constant.
    std::vector<std::vector<double>> assembleBoundaryMatrix(double c) {
        int nNodes = nodes.size();
        std::vector<std::vector<double>> B(nNodes, std::vector<double>(nNodes, 0.0));

        // Loop over each boundary edge.
        for (const auto& edge : boundaries) {
            int i = edge[0], j = edge[1];
            // Compute edge length between node i and node j.
            double dx = nodes[j].x - nodes[i].x;
            double dy = nodes[j].y - nodes[i].y;
            double edgeLength = std::sqrt(dx * dx + dy * dy);

            // Compute contribution using the factor (c*edgeLength/6).
            double factor = c * edgeLength / 6.0;
            B[i][i] += factor * 2;
            B[i][j] += factor * 1;
            B[j][i] += factor * 1;
            B[j][j] += factor * 2;
        }
        return B;
    }

    // Assemble the boundary load vector G.
    // c is a scaling constant and u0 is the boundary load value.
    std::vector<double> assembleBoundaryVector(double c, double u0) {
        int nNodes = nodes.size();
        std::vector<double> G(nNodes, 0.0);

        for (const auto& edge : boundaries) {
            int i = edge[0], j = edge[1];
            double dx = nodes[j].x - nodes[i].x;
            double dy = nodes[j].y - nodes[i].y;
            double edgeLength = std::sqrt(dx * dx + dy * dy);

            double contribution = c * u0 * edgeLength / 2.0;
            G[i] += contribution;
            G[j] += contribution;
        }
        return G;
    }

    // Assemble the stiffness matrix A for constant thermal conductivity 'a'
    // (corresponds to IntMatrix.m).
    std::vector<std::vector<double>> assembleStiffnessMatrix(double a) {
        int nNodes = nodes.size();
        std::vector<std::vector<double>> A(nNodes, std::vector<double>(nNodes, 0.0));

        // Loop over each element (triangle).
        for (const auto& elem : elements) {
            // Get the three node indices.
            int i = elem[0], j = elem[1], k = elem[2];

            // Get node coordinates.
            double x1 = nodes[i].x, y1 = nodes[i].y;
            double x2 = nodes[j].x, y2 = nodes[j].y;
            double x3 = nodes[k].x, y3 = nodes[k].y;

            // Compute triangle area using the formula.
            double area = 0.5 * std::abs(x1 * (y2 - y3) +
                                           x2 * (y3 - y1) +
                                           x3 * (y1 - y2));

            // Compute the gradients of the shape functions.
            // Note: For a linear triangular element the gradients are constant.
            double factor = 1.0 / (2.0 * area);
            // gradN for each local node (first row: derivative with respect to x, second row: derivative with respect to y).
            std::array<double, 3> dNdx = { (y2 - y3) * factor, (y3 - y1) * factor, (y1 - y2) * factor };
            std::array<double, 3> dNdy = { (x3 - x2) * factor, (x1 - x3) * factor, (x2 - x1) * factor };

            // Compute local stiffness matrix: localA(i,j) = a*(gradN_i dot gradN_j)*area.
            std::array<std::array<double, 3>, 3> localA = {0};
            for (int m = 0; m < 3; ++m) {
                for (int n = 0; n < 3; ++n) {
                    localA[m][n] = a * (dNdx[m] * dNdx[n] + dNdy[m] * dNdy[n]) * area;
                }
            }

            // Assemble local matrix into global stiffness matrix.
            std::array<int, 3> localNodes = { i, j, k };
            for (int m = 0; m < 3; ++m) {
                for (int n = 0; n < 3; ++n) {
                    A[localNodes[m]][localNodes[n]] += localA[m][n];
                }
            }
        }
        return A;
    }

    // Assemble stiffness matrix with variable conductivity at nodes 
    // Here, 'aCoeff' is a vector of nodal values.
    std::vector<std::vector<double>> assembleStiffnessMatrix2(const std::vector<double>& aCoeff) {
        int nNodes = nodes.size();
        std::vector<std::vector<double>> A(nNodes, std::vector<double>(nNodes, 0.0));

        for (const auto& elem : elements) {
            int i = elem[0], j = elem[1], k = elem[2];
            double x1 = nodes[i].x, y1 = nodes[i].y;
            double x2 = nodes[j].x, y2 = nodes[j].y;
            double x3 = nodes[k].x, y3 = nodes[k].y;

            double area = 0.5 * std::abs(x1 * (y2 - y3) +
                                           x2 * (y3 - y1) +
                                           x3 * (y1 - y2));

            double factor = 1.0 / (2.0 * area);
            std::array<double, 3> dNdx = { (y2 - y3) * factor, (y3 - y1) * factor, (y1 - y2) * factor };
            std::array<double, 3> dNdy = { (x3 - x2) * factor, (x1 - x3) * factor, (x2 - x1) * factor };

            // Average the conductivity over the element.
            double a_elem = (aCoeff[i] + aCoeff[j] + aCoeff[k]) / 3.0;

            std::array<std::array<double, 3>, 3> localA = {0};
            for (int m = 0; m < 3; ++m) {
                for (int n = 0; n < 3; ++n) {
                    localA[m][n] = a_elem * (dNdx[m] * dNdx[n] + dNdy[m] * dNdy[n]) * area;
                }
            }

            std::array<int, 3> localNodes = { i, j, k };
            for (int m = 0; m < 3; ++m) {
                for (int n = 0; n < 3; ++n) {
                    A[localNodes[m]][localNodes[n]] += localA[m][n];
                }
            }
        }
        return A;
    }

    // Assemble the stiffness matrix using quadrature 
    // Although similar to the variable coefficient method, here we explicitly use
    // the centroid for quadrature.
    std::vector<std::vector<double>> assembleStiffnessMatrixQuadrature(const std::vector<double>& aCoeff) {
        int nNodes = nodes.size();
        std::vector<std::vector<double>> A(nNodes, std::vector<double>(nNodes, 0.0));

        for (const auto& elem : elements) {
            int i = elem[0], j = elem[1], k = elem[2];
            double x1 = nodes[i].x, y1 = nodes[i].y;
            double x2 = nodes[j].x, y2 = nodes[j].y;
            double x3 = nodes[k].x, y3 = nodes[k].y;

            double area = 0.5 * std::abs(x1 * (y2 - y3) +
                                           x2 * (y3 - y1) +
                                           x3 * (y1 - y2));

            double factor = 1.0 / (2.0 * area);
            std::array<double, 3> dNdx = { (y2 - y3) * factor, (y3 - y1) * factor, (y1 - y2) * factor };
            std::array<double, 3> dNdy = { (x3 - x2) * factor, (x1 - x3) * factor, (x2 - x1) * factor };

            // Compute the centroid of the triangle.
            double x_quad = (x1 + x2 + x3) / 3.0;
            double y_quad = (y1 + y2 + y3) / 3.0;
            // In this example the shape functions at the centroid are all 1/3.
            double a_quad = (aCoeff[i] + aCoeff[j] + aCoeff[k]) / 3.0;

            std::array<std::array<double, 3>, 3> localA = {0};
            for (int m = 0; m < 3; ++m) {
                for (int n = 0; n < 3; ++n) {
                    localA[m][n] = a_quad * (dNdx[m] * dNdx[n] + dNdy[m] * dNdy[n]) * area;
                }
            }
            std::array<int, 3> localNodes = { i, j, k };
            for (int m = 0; m < 3; ++m) {
                for (int n = 0; n < 3; ++n) {
                    A[localNodes[m]][localNodes[n]] += localA[m][n];
                }
            }
        }
        return A;
    }

    // Assemble the load vector F using a quadrature rule 
    // f_handle is a function taking x and y as parameters and returning a double.
    std::vector<double> assembleLoadVectorQuadrature(std::function<double(double, double)> f_handle) {
        int nNodes = nodes.size();
        std::vector<double> F(nNodes, 0.0);

        for (const auto& elem : elements) {
            int i = elem[0], j = elem[1], k = elem[2];
            double x1 = nodes[i].x, y1 = nodes[i].y;
            double x2 = nodes[j].x, y2 = nodes[j].y;
            double x3 = nodes[k].x, y3 = nodes[k].y;

            double area = 0.5 * std::abs(x1 * (y2 - y3) +
                                           x2 * (y3 - y1) +
                                           x3 * (y1 - y2));

            // Compute midpoints of each edge.
            double q1x = (x1 + x2) / 2.0, q1y = (y1 + y2) / 2.0;
            double q2x = (x2 + x3) / 2.0, q2y = (y2 + y3) / 2.0;
            double q3x = (x3 + x1) / 2.0, q3y = (y3 + y1) / 2.0;

            // Evaluate the source function f at the quadrature points.
            double f1 = f_handle(q1x, q1y);
            double f2 = f_handle(q2x, q2y);
            double f3 = f_handle(q3x, q3y);

            // Local load vector using the quadrature rule.
            double factor = area / 3.0;
            F[i] += factor * f1;
            F[j] += factor * f2;
            F[k] += factor * f3;
        }
        return F;
    }
};

