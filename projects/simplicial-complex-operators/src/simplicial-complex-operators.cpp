// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"
#include <Eigen/Sparse>
#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

typedef Eigen::Triplet<double> T; // should the type be size_t instead of double?

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    // wg: i can't seem to find the definiiton of requireVertexIndices anywhere, i wonder how this works
    // oh these are all things in geometry central i think 
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input: 
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {    
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();
    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    std::vector<T> tripletList;
    tripletList.reserve(2*mesh->nEdges()); // two non-zero entries per row
    /*
        to construct A0 we need to iterate over all the edges, get the indices of its vertices
        and then construct a tripe T(i,j, v_ij) -> T(edge_num, vert_1_num, 1), T(edge_num, vert_2_num, 1)
    */
    size_t edge_idx;
    size_t first_vert_idx;
    size_t second_vert_idx;

    for (Edge e : mesh->edges())
    {
        // do for in loops generally give pointers or references back?
        edge_idx = geometry->edgeIndices[e];
        first_vert_idx = geometry->vertexIndices[e.firstVertex()];
        second_vert_idx = geometry->vertexIndices[e.secondVertex()]; 

        tripletList.push_back(T(edge_idx, first_vert_idx, 1));
        tripletList.push_back(T(edge_idx, second_vert_idx, 1));
    }

    SparseMatrix<size_t> mat(mesh->nEdges(), mesh->nVertices());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input: assuming that we're only concerned with simplectic surfaces
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    
    /*
        iterate over the faces. get some half-edge
        keep traversing using the half edge, at each iteration
            get the corresponding edge object, get its index and put it in the result matrix
        probs need to keep the first edge so that we can tell where we started
    */
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    std::vector<T> tripletList;
    tripletList.reserve(3*mesh->nFaces());

    size_t f_idx;
    size_t e_idx;
    Edge currentEdge;
    Halfedge first_halfEdge;
    Halfedge current_halfEdge;
    for (Face f : mesh->faces())
    {
        f_idx = geometry->faceIndices[f];
        first_halfEdge = f.halfedge();
        current_halfEdge = f.halfedge();
        do  {
            currentEdge = current_halfEdge.edge();
            e_idx = geometry->edgeIndices[currentEdge];
            tripletList.push_back(T(f_idx, e_idx, 1));

            current_halfEdge = current_halfEdge.next();
        } while (current_halfEdge != first_halfEdge);
    }

    SparseMatrix<size_t> mat(mesh->nFaces(), mesh->nEdges());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    return Vector<size_t>::Zero(1);
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    return Vector<size_t>::Zero(1);
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    return Vector<size_t>::Zero(1);
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    return false; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}