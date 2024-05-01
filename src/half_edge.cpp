#include <include/half_edge.hpp>
#include <Eigen/Dense>

/*
#################################################################################
#                       Vertex-related Helper Functions                         #
#################################################################################
*/

// Optinal TODO: Iterate through all neighbour vertices around the vertex
// Helpful when you implement the average degree computation of the mesh
std::vector<std::shared_ptr<Vertex>> Vertex::neighbor_vertices() {
    std::vector<std::shared_ptr<Vertex>> neighborhood;
    auto half = this->he;
    do {
        neighborhood.push_back(half->twin->vertex);
        half = half->twin->next;
    } while (half != this->he);
    return neighborhood;
}

void HalfEdge::remove() {
    this->exists = false;

    vertex->he = prev->twin;
    twin->vertex->update_to_vetex(vertex);

    if (next->next == prev) {
        next->edge->exists = false;
        prev->exists = false;
        next->exists = false;
        face->exists = false;

        prev->twin->twin = next->twin;
        next->twin->twin = prev->twin;

        prev->edge->he = prev->twin;
        next->twin->edge = prev->edge;

        prev->vertex->he = next->twin;
    } else {
        prev->next = next;
        next->prev = prev;
    }

    twin->exists = false;
    if (twin->next->next == twin->prev) {

        twin->prev->edge->exists = false;
        twin->prev->exists = false;
        twin->next->exists = false;
        twin->face->exists = false;

        twin->prev->twin->twin = twin->next->twin;
        twin->next->twin->twin = twin->prev->twin;

        twin->next->edge->he = twin->next->twin;
        twin->prev->twin->edge = twin->next->edge;

        twin->prev->vertex->he = twin->next->twin;

    } else {
        twin->prev->next = twin->next;
        twin->next->prev = twin->prev;
    }
}


// TODO: Iterate through all half edges pointing away from the vertex
std::vector<std::shared_ptr<HalfEdge>> Vertex::neighbor_half_edges() {
    std::vector<std::shared_ptr<HalfEdge>> neighborhood;
    auto half = this->he;

    do {
        neighborhood.push_back(half);
        half = half->twin->next;
    } while (half != this->he);

    return neighborhood;
}


// TODO: Computate quadratic error metrics coefficient, which is a 5-d vector associated with each vertex
/*
    HINT:
        Please refer to homework description about how to calculate each element in the vector.
        The final results is stored in class variable "this->qem_coff"
*/
void Vertex::compute_qem_coeff() {
    this->qem_coff = Eigen::VectorXf(5);
    auto neighbors = this->neighbor_vertices();
    float q1 = 0, q5 = 0; // this is stupid;
    q1 = neighbors.size();
    Eigen::Vector3f vertex_pos_sum = {0, 0, 0};

    for (const auto vertex: neighbors) {
        assert(vertex != nullptr);
        vertex_pos_sum += vertex->pos;
        q5 += vertex->pos.dot(vertex->pos);
    }

    this->qem_coff << q1, vertex_pos_sum.x(), vertex_pos_sum.y(), vertex_pos_sum.z(), q5;
}

void Vertex::update_to_vetex(std::shared_ptr<Vertex> new_vertex) {
    this->exists = false;
    this->he->vertex = new_vertex;

    for (auto half_edge: this->neighbor_half_edges()) {
        half_edge->vertex = new_vertex;
    }
}

/*
#################################################################################
#                         Face-related Helper Functions                         #
#################################################################################
*/

// TODO: Iterate through all member vertices of the face
std::vector<std::shared_ptr<Vertex>> Face::vertices() {
    std::vector<std::shared_ptr<Vertex>> face_vertices;
    auto half = this->he;
    do {
        face_vertices.push_back(half->vertex);
        half = half->next;
    } while (half != this->he);
    return face_vertices;
}


// TODO: implement this function to compute the area of the triangular face
float Face::get_area() {
    auto vetices = this->vertices();
    // calcaulte the area of the triangle
    auto p1 = vetices[0]->pos, p2 = vetices[1]->pos, p3 = vetices[2]->pos;
    auto v1 = p2 - p1, v2 = p3 - p1;
    auto cross = v1.cross(v2);
    return 0.5 * cross.norm(); // 0.5 * |v1 x v2|
}

// TODO: implement this function to compute the signed volume of the triangular face
// reference: http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf eq.(5)
float Face::get_signed_volume() {
    auto vertices = this->vertices();
    auto p1 = vertices[0]->pos, p2 = vertices[1]->pos, p3 = vertices[2]->pos;
    return p1.dot(p2.cross(p3)) / 6;
}


/*
#################################################################################
#                         Edge-related Helper Functions                         #
#################################################################################
*/

/*
    TODO: 
        Compute the contraction information for the edge (v1, v2), which will be used later to perform edge collapse
            (i) The optimal contraction target v*
            (ii) The quadratic error metrics QEM, which will become the cost of contracting this edge
        The final results is stored in class variable "this->verts_contract_pos" and "this->qem"
    Please refer to homework description for more details
*/
void Edge::compute_contraction() {
    // compute the optimal contraction target v* and the QEM
    auto v1 = this->he->vertex, v2 = this->he->twin->vertex;
    auto q1 = v1->qem_coff, q2 = v2->qem_coff;
    auto q = q1 + q2;

    const auto calc_new_pos = [&, this]() {
        auto total_pos_sum = q1.segment(1, 3) + q2.segment(1, 3);
        auto weight = q1(0) + q2(0);
        return total_pos_sum / weight;
    };

    const auto calc_error = [](const Eigen::Vector3f &pos, const Eigen::VectorXf &q) {
        return q(0) * pos.dot(pos) - 2 * pos.dot(q.segment(1, 3)) + q(4);
    };

    this->verts_contract_pos = calc_new_pos();
    this->qem = calc_error(this->verts_contract_pos, q);
}

/*
    TODO: 
        Perform edge contraction functionality, which we write as (v1 ,v2) → v*, 
            (i) Moves the vertex v1 to the new position v*, remember to update all corresponding attributes,
            (ii) Connects all incident edges of v1 and v2 to v*, and remove the vertex v2,
            (iii) All faces, half edges, and edges associated with this collapse edge will be removed.
    HINT: 
        (i) Pointer reassignments
        (ii) When you want to remove mesh components, simply set their "exists" attribute to False
    Please refer to homework description for more details
*/
void Edge::edge_contraction() {
    // move the vertex v1 to the new position v*
    this->exists = false;

    auto v1 = this->he->vertex, v2 = this->he->twin->vertex;
    auto combine_q = v1->qem_coff + v2->qem_coff;

    // this operation will remove v2 from the mesh
    he->remove();

    v1->pos = this->verts_contract_pos;
    v1->qem_coff = combine_q;
}