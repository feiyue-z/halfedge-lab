#include "MeshDataStructure.h"

#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <unordered_set>

using namespace Eigen;
using namespace std;

void MeshDataStructure::init(const Mesh &mesh)
{
    map<PointPair, Halfedge*, PointPairComparator> halfedgeMap;

    // populate vertices
    for (const Vector3f &v : mesh._vertices) {
        Vertex *vertex = new Vertex(v);
        _vertices.push_back(vertex);
    }

    // assume triangle face
    for (const Vector3i &f : mesh._faces) {
        // face
        Face *face = new Face;
        _faces.push_back(face);

        Halfedge* faceHalfedges[3]; // halfedges in this face

        for (int i = 0; i < 3; i++) {
            // halfedge
            Halfedge *he = new Halfedge;
            _halfedges.push_back(he);
            faceHalfedges[i] = he;

            // next
            if (i != 0) {
                faceHalfedges[i - 1]->next = he;
            }
            if (i == 2) {
                he->next = faceHalfedges[0];
            }

            // twin
            PointPair pp(mesh._vertices[f[i]], mesh._vertices[f[(i + 1) % 3]]);
            if (halfedgeMap.contains(pp)) {
                he->twin = halfedgeMap.at(pp);
                he->twin->twin = he;
            } else {
                halfedgeMap[pp] = he;
            }

            // vertex
            he->vertex = _vertices[f[i]];
            if (he->vertex->halfedge == nullptr) {
                he->vertex->halfedge = he;
            }

            // edge
            if (he->twin == nullptr) {
                Edge *edge = new Edge;
                _edges.push_back(edge);
                edge->halfedge = he;
                he->edge = edge;
            } else {
                he->edge = he->twin->edge;
            }

            // face
            he->face = face;
        }

        face->halfedge = faceHalfedges[0];
    }
}

void MeshDataStructure::validate()
{
    for (Halfedge *he : _halfedges) {
        // each halfedge is not missing data
        assert(he->twin != nullptr);
        assert(he->next != nullptr);
        assert(he->vertex != nullptr);
        assert(he->edge != nullptr);
        assert(he->face != nullptr);

        // each halfedge's twin's twin points back to itself
        assert(he->twin->twin == he);

        // each halfedge's vertex and that of its twin is different
        assert(he->vertex != he->twin->vertex);

        // traverse each halfedge's next until it points back to itself
        unordered_set<Halfedge*> visited;
        Halfedge *curr = he;
        do {
            // detect loop in traversal
            assert(visited.find(curr) == visited.end());

            visited.insert(curr);
            curr = curr->next;
        }
        while (curr != he);
        assert(curr == he);
    }

    for (Vertex *v : _vertices) {
        // each vertex is not missing data
        assert(v->halfedge != nullptr);

        // each vertex's halfedge's vertex points back to itself
        assert(v->halfedge->vertex == v);

        // traverse halfedges around vertex until it points back to itself
        unordered_set<Halfedge*> visited;
        Halfedge *curr = v->halfedge;
        do {
            // detect loop in traversal
            assert(visited.find(curr) == visited.end());

            visited.insert(curr);
            curr = curr->twin->next;
        }
        while (curr != v->halfedge);
        assert(curr == v->halfedge);
    }

    for (Edge *e : _edges) {
        // each edge is not missing data
        assert(e->halfedge != nullptr);
    }

    for (Face *f : _faces) {
        // each face is not missing data
        assert(f->halfedge != nullptr);
    }
}

// n: vertex degree
float getU(int n)
{
    if (n == 3) {
        return 3.f / 16.f;
    }
    return (5.f / 8.f - pow((3.f / 8.f + (cos(2.f * M_PI / n)) / 4.f), 2)) / n;
}

void MeshDataStructure::subdivide(int numIterations)
{
    for (int i = 0; i < numIterations; i++) {
        cout << "i = " << i << endl;
        subdivideOnce();
        validate();
    }
}

void MeshDataStructure::subdivideOnce()
{
    vector<Halfedge*> newHalfedges;
    vector<Edge*> newEdges;
    vector<Face*> newFaces;

    map<Edge*, Vertex*> verticesSplitFromEdge;
    map<PointPair, Halfedge*, PointPairComparator> halfedgeMap;
    map<Vertex*, Vector3f> updatedPoints;

    // calculate new points for existing vertices
    for (Vertex *v : _vertices) {
        int cnt = 0;
        Vector3f sum(0, 0, 0);
        Halfedge *curr = v->halfedge;

        do {
            sum += curr->next->vertex->point;
            curr = curr->twin->next;
            cnt++;
        }
        while (curr != v->halfedge);

        float u = getU(cnt);
        Vector3f point = (1 - cnt * u) * v->point + u * sum;
        updatedPoints[v] = point;
    }

    for (Face *f : _faces) {
        // halfedges in this face
        // in CCW order
        Halfedge* he1 = f->halfedge;
        Halfedge* he2 = f->halfedge->next;
        Halfedge* he3 = f->halfedge->next->next;

        // split edge if this edge has not been split
        if (!verticesSplitFromEdge.contains(he1->edge)) {
            splitEdge(he1->edge, verticesSplitFromEdge);
        }
        if (!verticesSplitFromEdge.contains(he2->edge)) {
            splitEdge(he2->edge, verticesSplitFromEdge);
        }
        if (!verticesSplitFromEdge.contains(he3->edge)) {
            splitEdge(he3->edge, verticesSplitFromEdge);
        }

        // new vertices from splitting the edges in this face
        Vertex *nv1 = verticesSplitFromEdge[he1->edge];
        Vertex *nv2 = verticesSplitFromEdge[he2->edge];
        Vertex *nv3 = verticesSplitFromEdge[he3->edge];

        // sub-triangle 1
        vector<Vertex *> faceVertices = {nv1, nv2, nv3};
        updateStructure(faceVertices,
                        newHalfedges,
                        newEdges,
                        newFaces,
                        halfedgeMap);

        // sub-triangle 2
        Vertex *v2 = f->halfedge->next->vertex;
        faceVertices = {nv1, v2, nv2};
        updateStructure(faceVertices,
                        newHalfedges,
                        newEdges,
                        newFaces,
                        halfedgeMap);

        // sub-triangle 3
        Vertex *v3 = f->halfedge->next->next->vertex;
        faceVertices = {nv2, v3, nv3};
        updateStructure(faceVertices,
                        newHalfedges,
                        newEdges,
                        newFaces,
                        halfedgeMap);

        // sub-triangle 4
        Vertex *v1 = f->halfedge->vertex;
        faceVertices = {nv3, v1, nv1};
        updateStructure(faceVertices,
                        newHalfedges,
                        newEdges,
                        newFaces,
                        halfedgeMap);
    }

    // update existing vertices with pre-calculated points
    for (Vertex *v : _vertices) {
        v->point = updatedPoints[v];
    }

    // add new vertices
    for (const auto &pair : verticesSplitFromEdge) {
        _vertices.push_back(pair.second);
    }

    // replace old structures with new ones
    _edges = newEdges;
    _faces = newFaces;
    _halfedges = newHalfedges;
}

void MeshDataStructure::splitEdge(Edge *e, map<Edge*, Vertex*> &verticesSplitFromEdge)
{
    // adjacent vertices
    Vertex *adj1 = e->halfedge->vertex;
    Vertex *adj2 = e->halfedge->twin->vertex;

    // opposite vertices
    Vertex *opp1 = e->halfedge->next->next->vertex;
    Vertex *opp2 = e->halfedge->twin->next->next->vertex;

    // new vertex
    Vector3f point = 3.f * adj1->point / 8.f +
                     3.f * adj2->point / 8.f +
                     opp1->point / 8.f +
                     opp2->point / 8.f;
    Vertex *v = new Vertex(point);
    verticesSplitFromEdge[e] = v;
}

// 3 vertices of a face
// generates new halfedge, edge, face in the given face
void MeshDataStructure::updateStructure(const vector<Vertex*> &faceVertices,
                                        vector<Halfedge*> &halfedges,
                                        vector<Edge*> &edges,
                                        vector<Face*> &faces,
                                        map<PointPair, Halfedge*, PointPairComparator> &halfedgeMap)
{
    Face *face = new Face;
    faces.push_back(face);

    Halfedge* faceHalfedges[3]; // halfedges in this face

    for (int i = 0; i < 3; i++) {
        // halfedge
        Halfedge *he = new Halfedge;
        halfedges.push_back(he);
        faceHalfedges[i] = he;

        // next
        if (i != 0) {
            faceHalfedges[i - 1]->next = he;
        }
        if (i == 2) {
            he->next = faceHalfedges[0];
        }

        // twin
        PointPair pp(faceVertices[i]->point, faceVertices[(i + 1) % 3]->point);
        if (halfedgeMap.contains(pp)) {
            he->twin = halfedgeMap.at(pp);
            he->twin->twin = he;
        } else {
            halfedgeMap[pp] = he;
        }

        // vertex
        he->vertex = faceVertices[i];
        he->vertex->halfedge = he;

        // edge
        if (he->twin == nullptr) {
            Edge *edge = new Edge;
            edges.push_back(edge);
            edge->halfedge = he;
            he->edge = edge;
        } else {
            he->edge = he->twin->edge;
        }

        // face
        he->face = face;
    }

    face->halfedge = faceHalfedges[0];
}

void MeshDataStructure::simplify(int numFaceRemoves)
{
    map<Vertex*, Matrix4f> vertexQuadricMap;
    for (const Face *f : _faces) {
        computeVertexQuadric(f, vertexQuadricMap);
    }

    multiset<Edge*, EdgeCostComparator> edgeCostSet;
    for (Edge *e : _edges) {
        e->cost = getEdgeCostAndOptimalVertex(e, vertexQuadricMap);
        edgeCostSet.insert(e);
    }

    int numEdgeRemoves = round(numFaceRemoves / 2);
    for (int i = 0; i < numEdgeRemoves; i++) {
        auto it = edgeCostSet.begin();
        Edge *eLowestCost = *it;

        eraseEdgeFromSet(eLowestCost, edgeCostSet);
        if (!canCollapse(eLowestCost)) {
            i--;
            continue;
        }

        collapseAndUpdateCost(eLowestCost, vertexQuadricMap, edgeCostSet);
        validate();
    }
}

void MeshDataStructure::eraseEdgeFromSet(Edge *e, multiset<Edge*, EdgeCostComparator> &edgeCostSet)
{
    for (auto it = edgeCostSet.begin(); it != edgeCostSet.end(); it++) {
        if (*it == e) {
            edgeCostSet.erase(it);
            return;
        }
    }
    cout << "not found" << endl;
}

bool MeshDataStructure::edgeIsInSet(Edge *e, multiset<Edge*, EdgeCostComparator> &edgeCostSet)
{
    for (auto it = edgeCostSet.begin(); it != edgeCostSet.end(); it++) {
        if (*it == e) {
            return true;
        }
    }
    return false;
}

bool MeshDataStructure::canCollapse(Edge *e)
{
    // e might be deleted during previous collapse operations
    // ???
    auto it = find(_edges.begin(), _edges.end(), e);
    if (it == _edges.end()) {
        // not found in _edges
        return false;
    }

    Vertex *v1 = e->halfedge->vertex;
    Vertex *v2 = e->halfedge->twin->vertex;

    unordered_set<Vertex*> neighbors1 = getNeighborVertices(v1);
    unordered_set<Vertex*> neighbors2 = getNeighborVertices(v2);
    unordered_set<Vertex*> shared;

    if (neighbors1.size() < neighbors2.size()) {
        for (const auto& each : neighbors1) {
            if (neighbors2.find(each) != neighbors2.end()) {
                shared.insert(each);
            }
        }
    } else {
        for (const auto& each : neighbors2) {
            if (neighbors1.find(each) != neighbors1.end()) {
                shared.insert(each);
            }
        }
    }

    if (shared.size() > 2) {
        return false;
    }

    for (const auto& each : shared) {
        unordered_set<Vertex*> neighbors = getNeighborVertices(each);

        if (neighbors.size() == 3) {
            return false;
        }
    }
    return true;
}

unordered_set<MeshDataStructure::Vertex*> MeshDataStructure::getNeighborVertices(const Vertex *v)
{
    unordered_set<Vertex*> neighbors;
    Halfedge *curr = v->halfedge;

    do {
        neighbors.insert(curr->twin->vertex);
        curr = curr->twin->next;
    }
    while (curr != v->halfedge);

    return neighbors;
}

void MeshDataStructure::computeVertexQuadric(const Face *f,
                                             map<Vertex*, Matrix4f> &vertexQuadricMap)
{
    Vertex *v1 = f->halfedge->vertex;
    Vertex *v2 = f->halfedge->next->vertex;
    Vertex *v3 = f->halfedge->next->next->vertex;

    // normal of this face
    Vector3f N = getFaceNormal(v1, v2, v3);

    // plane offset
    float d = -N.dot(v1->point);

    Vector4f v(N.x(), N.y(), N.z(), d);

    // plane quadric
    Matrix4f Q = v * v.transpose();

    // sum quadrics
    vertexQuadricMap[v1] += Q;
    vertexQuadricMap[v2] += Q;
    vertexQuadricMap[v3] += Q;
}

Vector3f MeshDataStructure::getFaceNormal(Vertex *v1, Vertex *v2, Vertex *v3)
{
    return (v2->point - v1->point).cross(v3->point - v1->point).normalized();
}

float MeshDataStructure::getEdgeCostAndOptimalVertex(Edge *e,
                                                     map<Vertex*, Matrix4f> &vertexQuadricMap)
{
    // 2 vertices on this edge
    Vertex *v1 = e->halfedge->vertex;
    Vertex *v2 = e->halfedge->twin->vertex;

    // quadric of optical vertex
    Matrix4f Q = vertexQuadricMap.at(v1) + vertexQuadricMap.at(v2);

    // solve for optical vertex
    Matrix3f Q_prime = Q.block<3, 3>(0, 0);
    Vector3f b = -Q.block<3, 1>(0, 3);
    Vector3f v_prime = Q_prime.ldlt().solve(b);

    Vector4f v_optimal(v_prime[0], v_prime[1], v_prime[2], 1.f);
    e->optimalVertex = new Vertex(Vector3f(v_prime));
    vertexQuadricMap[e->optimalVertex] = Q;

    // edge cost
    float cost = v_optimal.transpose() * Q * v_optimal;
    return cost;
}

void MeshDataStructure::collapse(Edge *eCollap)
{
    Halfedge *he1a = eCollap->halfedge;
    Halfedge *he1b = eCollap->halfedge->next;
    Halfedge *he1c = eCollap->halfedge->next->next;

    Halfedge *he2a = eCollap->halfedge->twin;
    Halfedge *he2b = eCollap->halfedge->twin->next;
    Halfedge *he2c = eCollap->halfedge->twin->next->next;

    // 2 vertices on this edge
    Vertex *v1 = he1a->vertex;
    Vertex *v2 = he2a->vertex;

    // 3 degenerate edges to remove
    unordered_set<Edge*> eRemove;
    eRemove.insert(eCollap);
    eRemove.insert(he1c->edge);
    eRemove.insert(he2c->edge);

    // 2 degenerate faces to remove
    Face *fToRemove1 = he1a->face;
    Face *fToRemove2 = he2a->face;

    // 6 degenerate halfedges to remove
    unordered_set<Halfedge*> heRemove1 = getHalfedgesInTriangle(he1a);
    unordered_set<Halfedge*> heRemove2 = getHalfedgesInTriangle(he2a);

    // new vertex in the middle of v1 and v2
    Vertex *vMid = new Vertex((v1->point + v2->point) / 2);
    _vertices.push_back(vMid);

    // update collapsed vertices with new vertex
//    unordered_set<Edge*> eRecalc1 = updateCollapsedVertex(eCollap, v1, vMid);
//    unordered_set<Edge*> eRecalc2 = updateCollapsedVertex(eCollap, v2, vMid);
    updateCollapsedVertex(eCollap, v1, vMid);
    updateCollapsedVertex(eCollap, v2, vMid);

    ensureValidHalfedgeForVertex(he1c->vertex, heRemove1);
    ensureValidHalfedgeForVertex(he2c->vertex, heRemove2);

    mergeCollapsedEdges(he1a);
    mergeCollapsedEdges(he2a);

    vMid->halfedge = he1b->twin->twin;

    // remove collapsed vertices
    int size1 = _vertices.size();
    _vertices.erase(remove(_vertices.begin(), _vertices.end(), v1), _vertices.end());
    _vertices.erase(remove(_vertices.begin(), _vertices.end(), v2), _vertices.end());
    int size2 = _vertices.size();
    assert(size1 - size2 == 2);

    // remove collapsed edge
    size1 = _edges.size();
    for (Edge *each : eRemove) {
        _edges.erase(remove(_edges.begin(), _edges.end(), each), _edges.end());
    }
    size2 = _edges.size();
    assert(size1 - size2 == 3);

    // remove collapsed faces
    size1 = _faces.size();
    _faces.erase(remove(_faces.begin(), _faces.end(), fToRemove1), _faces.end());
    _faces.erase(remove(_faces.begin(), _faces.end(), fToRemove2), _faces.end());
    size2 = _faces.size();
    assert(size1 - size2 == 2);

    // remove collapsed halfedges
    size1 = _halfedges.size();
    for (Halfedge* each : heRemove1) {
        _halfedges.erase(remove(_halfedges.begin(), _halfedges.end(), each), _halfedges.end());
    }
    for (Halfedge* each : heRemove2) {
        _halfedges.erase(remove(_halfedges.begin(), _halfedges.end(), each), _halfedges.end());
    }
    size2 = _halfedges.size();
    assert(size1 - size2 == 6);

    delete v1;
    delete v2;
    delete fToRemove1;
    delete fToRemove2;
    for (Edge *each : eRemove) {
        delete each;
    }
    for (Halfedge* each : heRemove1) {
        delete each;
    }
    for (Halfedge* each : heRemove2) {
        delete each;
    }
}

void MeshDataStructure::collapseAndUpdateCost(Edge *eCollap,
                                              map<Vertex*, Matrix4f> &vertexQuadricMap,
                                              multiset<Edge*, EdgeCostComparator> &edgeCostSet)
{
    Halfedge *he1a = eCollap->halfedge;
    Halfedge *he1b = eCollap->halfedge->next;
    Halfedge *he1c = eCollap->halfedge->next->next;

    Halfedge *he2a = eCollap->halfedge->twin;
    Halfedge *he2b = eCollap->halfedge->twin->next;
    Halfedge *he2c = eCollap->halfedge->twin->next->next;

    // 2 vertices on this edge
    Vertex *v1 = he1a->vertex;
    Vertex *v2 = he2a->vertex;

    // 3 degenerate edges to remove
    unordered_set<Edge*> eRemove;
    eRemove.insert(eCollap);
    eRemove.insert(he1c->edge);
    eRemove.insert(he2c->edge);

    // 2 degenerate faces to remove
    Face *fToRemove1 = he1a->face;
    Face *fToRemove2 = he2a->face;

    // 6 degenerate halfedges to remove
    unordered_set<Halfedge*> heRemove1 = getHalfedgesInTriangle(he1a);
    unordered_set<Halfedge*> heRemove2 = getHalfedgesInTriangle(he2a);

    // new vertex in the middle of v1 and v2
    Vertex *vMid = eCollap->optimalVertex;
    _vertices.push_back(vMid);

    // update collapsed vertices with new vertex
    unordered_set<Edge*> eRecalc1 = updateCollapsedVertex(eCollap, v1, vMid);
    unordered_set<Edge*> eRecalc2 = updateCollapsedVertex(eCollap, v2, vMid);

    // calculate new cost for affected edges
    for (Edge *e : eRecalc1) {
        e->cost = getEdgeCostAndOptimalVertex(e, vertexQuadricMap);
        eraseEdgeFromSet(e, edgeCostSet);
        edgeCostSet.insert(e);
    }
    for (Edge *e : eRecalc2) {
        e->cost = getEdgeCostAndOptimalVertex(e, vertexQuadricMap);
        eraseEdgeFromSet(e, edgeCostSet);
        edgeCostSet.insert(e);
    }

    ensureValidHalfedgeForVertex(he1c->vertex, heRemove1);
    ensureValidHalfedgeForVertex(he2c->vertex, heRemove2);

    mergeCollapsedEdges(he1a);
    mergeCollapsedEdges(he2a);

    vMid->halfedge = he1b->twin->twin;

    // remove collapsed vertices
    int size1 = _vertices.size();
    _vertices.erase(remove(_vertices.begin(), _vertices.end(), v1), _vertices.end());
    _vertices.erase(remove(_vertices.begin(), _vertices.end(), v2), _vertices.end());
    int size2 = _vertices.size();
    assert(size1 - size2 == 2);

    // remove errors of collapsed vertices
    size1 = vertexQuadricMap.size();
    vertexQuadricMap.erase(v1);
    vertexQuadricMap.erase(v2);
    size2 = vertexQuadricMap.size();
    assert(size1 - size2 == 2);

    // remove collapsed edge
    size1 = _edges.size();
    for (Edge *each : eRemove) {
        _edges.erase(remove(_edges.begin(), _edges.end(), each), _edges.end());
    }
    size2 = _edges.size();
    assert(size1 - size2 == 3);

    // remove collapsed faces
    size1 = _faces.size();
    _faces.erase(remove(_faces.begin(), _faces.end(), fToRemove1), _faces.end());
    _faces.erase(remove(_faces.begin(), _faces.end(), fToRemove2), _faces.end());
    size2 = _faces.size();
    assert(size1 - size2 == 2);

    // remove collapsed halfedges
    size1 = _halfedges.size();
    for (Halfedge* each : heRemove1) {
        _halfedges.erase(remove(_halfedges.begin(), _halfedges.end(), each), _halfedges.end());
    }
    for (Halfedge* each : heRemove2) {
        _halfedges.erase(remove(_halfedges.begin(), _halfedges.end(), each), _halfedges.end());
    }
    size2 = _halfedges.size();
    assert(size1 - size2 == 6);

    delete v1;
    delete v2;
    delete fToRemove1;
    delete fToRemove2;
    for (Edge *each : eRemove) {
        delete each;
    }
    for (Halfedge* each : heRemove1) {
        delete each;
    }
    for (Halfedge* each : heRemove2) {
        delete each;
    }
}

// keeps edge a, remove edge b
void MeshDataStructure::mergeCollapsedEdges(const Halfedge *he)
{
    Halfedge *aIn = he->next;
    Halfedge *bIn = he->next->next;

    Halfedge *aOut = aIn->twin;
    Halfedge *bOut = bIn->twin;

    bOut->twin = aOut;
    aOut->twin = bOut;

    aOut->edge->halfedge = aOut;
    bOut->edge->halfedge = nullptr;
    bOut->edge = aOut->edge;
}

void MeshDataStructure::ensureValidHalfedgeForVertex(Vertex *v, unordered_set<Halfedge*> heRemove)
{
    auto it = heRemove.find(v->halfedge);
    if (it != heRemove.end()) {
        v->halfedge = v->halfedge->twin->next;
    }
}

unordered_set<MeshDataStructure::Halfedge*> MeshDataStructure::getHalfedgesInTriangle(Halfedge *he)
{
    unordered_set<Halfedge*> halfedges;
    Halfedge* curr = he;

    do {
        halfedges.insert(curr);
        curr = curr->next;
    }
    while (curr != he);

    return halfedges;
}

unordered_set<MeshDataStructure::Edge*> MeshDataStructure::updateCollapsedVertex(const Edge* e,
                                                                                 const Vertex* vOld,
                                                                                 Vertex* vNew)
{
    unordered_set<Edge*> affectedEdges;
    Halfedge *curr = vOld->halfedge;

    do {
        if (curr->edge != e) {
            curr->vertex = vNew;
            affectedEdges.insert(curr->edge);

            if (vNew->halfedge == nullptr) {
                vNew->halfedge = curr;
            }
        }

        curr = curr->twin->next;
    }
    while (curr != vOld->halfedge);

    return affectedEdges;
}

void MeshDataStructure::noise(int numIterations, float scale)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(-1, 1);

    for (int i = 0; i < numIterations; i++) {
        for (Vertex *v : _vertices) {
            v->point = v->point + Vector3f(1, 1, 1) * dist(gen) * scale;
        }
    }
}

void MeshDataStructure::denoise(int numIterations, float sigma_c, float sigma_s, float rho)
{
    for (int i = 0; i < numIterations; i++) {
        map<Vertex*, Vector3f> denoisedVertexValues;

        for (Vertex *v : _vertices) {
            Vector3f point = getDenoisedPoint(v, sigma_c, sigma_s, rho);
            denoisedVertexValues[v] = point;
        }

        for (Vertex *v : _vertices) {
            v->point = denoisedVertexValues.at(v);
        }
    }
    validate();
}

Vector3f MeshDataStructure::getDenoisedPoint(const Vertex* v,
                                         float sigma_c,
                                         float sigma_s,
                                         float rho)
{
    Vector3f N = getVertexNormal(v);
    unordered_set<Vertex*> neighbors = getNeighborVertices(v);

    float sum = 0.f;
    float normalizer = 0.f;

    for (Vertex *q : neighbors) {
        // distance
        float t = (v->point - q->point).norm();
        cout << "t = " << t << endl;
        if (t > rho) {
            continue;
        }

        //
        float h = N.dot(q->point - v->point);
        float w_c = exp(-t * t / (2 * sigma_c * sigma_c));
        float w_s = exp(-h * h / (2 * sigma_s * sigma_s));

        sum += w_c * w_s * h;
        normalizer += w_c * w_s;
    }

    if (normalizer != 0) {
        Vector3f x = N * (sum / normalizer);
        return v->point + N * (sum / normalizer);
    }
    return v->point;
}

void MeshDataStructure::remesh(int numIterations, float weight)
{
    float L = getMeanEdgeLength();

    cout << "4 * L / 3 = " << 4 * L / 3 << endl;
    cout << "4 * L / 5 = " << 4 * L / 5 << endl;

    // split
//    vector<Edge*> edgesToSplit;
//    for (Edge* e : _edges) {
//        float len = getEdgeLength(e);
//        if (len > 4 * L / 3) {
//            edgesToSplit.push_back(e);
//        }
//    }

//    map<Edge*, Vertex*> verticesSplitFromEdge;
//    for (Edge *e : edgesToSplit) {
//        cout << "split e = " << e << endl;
//        splitEdge(e, verticesSplitFromEdge);
//        splitTwoTrianglesIntoFour(e, verticesSplitFromEdge.at(e));
//    }

    // collapse
//    vector<Edge*> edgesToCollapse;
//    for (Edge *e : _edges) {
//        float len = getEdgeLength(e);
//        if (len < 4 * L / 5 && canCollapse(e)) {
//            edgesToCollapse.push_back(e);
//        }
//    }

//    int i = 0;
//    for (Edge *e : edgesToCollapse) {
//        // check because the pre-calculated edge-to-collapse
//        // may have been collapsed in previous operation
//        auto it = find(_edges.begin(), _edges.end(), e);
//        if (it != _edges.end()) {
//            cout << "collapse e = " << e << endl;
//            collapse(e);
//        }
//    }

    // flip
    for (Edge *e : _edges) {
        if (validToFlipForRemesh(e) && canFlip(e)) {
            flip(e);
        }
    }

    // calculate centroid
    map<Vertex*, Vector3f> centroidMap;
    for (Vertex *v : _vertices) {
        centroidMap[v] = getCentroid(v);
    }

    // move vertices tangentially to centroid
    for (Vertex *v : _vertices) {
        moveVertexTangentially(v, centroidMap.at(v), weight);
    }
}

float MeshDataStructure::getMeanEdgeLength()
{
    float sum = 0.f;
    for (Edge *e : _edges) {
        sum += getEdgeLength(e);
    }
    return sum / _edges.size();
}

float MeshDataStructure::getEdgeLength(Edge *e)
{
    Vector3f p1 = e->halfedge->vertex->point;
    Vector3f p2 = e->halfedge->twin->vertex->point;
    return (p1 - p2).norm();
}

void MeshDataStructure::splitTwoTrianglesIntoFour(Edge *e, Vertex *vMid)
{
    // halfedges in triangle 1
    Halfedge *he1 = e->halfedge;
    Halfedge *he2 = he1->next;
    Halfedge *he3 = he2->next;

    // halfedges in triangle 2
    Halfedge *he4 = e->halfedge->twin;
    Halfedge *he5 = he4->next;
    Halfedge *he6 = he5->next;

    Vertex *va = he1->vertex;
    Vertex *vb = he2->vertex;
    Vertex *vc = he3->vertex;
    Vertex *vd = he6->vertex;

    Face *face1 = he1->face;
    Face *face2 = he4->face;

    // split from face1
    Face *newFace1 = new Face;
    // split from face2
    Face *newFace2 = new Face;

    // split from e
    Edge *newE1 = new Edge;
    // new edges in another direction
    Edge *newE2 = new Edge;
    Edge *newE3 = new Edge;

    // split from he1
    Halfedge *newHe1 = new Halfedge;
    // split from he4
    Halfedge *newHe2 = new Halfedge;
    // new halfedges in another direction
    Halfedge *newHe3 = new Halfedge;
    Halfedge *newHe4 = new Halfedge;
    Halfedge *newHe5 = new Halfedge;
    Halfedge *newHe6 = new Halfedge;

    // update edge
    e->halfedge = he1;
//    cout << "e = " << e << endl;

    // populate edge
    newE1->halfedge = newHe1;
    newE2->halfedge = newHe3;
    newE3->halfedge = newHe4;

    // update halfedge
    he1->twin = newHe2;
    he4->twin = newHe1;

    he1->next = newHe5;
    he2->next = newHe3;
    he4->next = newHe6;
    he5->next = newHe4;

    he2->face = newFace1;
    he5->face = newFace2;

    // populate halfedge
    newHe1->twin = he4;
    newHe1->next = he2;
    newHe1->vertex = vMid;
    newHe1->edge = newE1;
    newHe1->face = newFace1;

    newHe2->twin = he1;
    newHe2->next = he5;
    newHe2->vertex = vMid;
    newHe2->edge = e;
    newHe2->face = newFace2;

    newHe3->twin = newHe5;
    newHe3->next = newHe1;
    newHe3->vertex = vc;
    newHe3->edge = newE2;
    newHe3->face = newFace1;

    newHe4->twin = newHe6;
    newHe4->next = newHe2;
    newHe4->vertex = vd;
    newHe4->edge = newE3;
    newHe4->face = newFace2;

    newHe5->twin = newHe3;
    newHe5->next = he3;
    newHe5->vertex = vMid;
    newHe5->edge = newE2;
    newHe5->face = face1;

    newHe6->twin = newHe4;
    newHe6->next = he6;
    newHe6->vertex = vMid;
    newHe6->edge = newE3;
    newHe6->face = face2;

    // populate face
    newFace1->halfedge = newHe1;
    newFace2->halfedge = newHe2;

    // populate vertex
    vMid->halfedge = newHe1;

    // add new data
    _vertices.push_back(vMid);

    _edges.push_back(newE1);
    _edges.push_back(newE2);
    _edges.push_back(newE3);

    // add to temporary container instead of directly modifying _edges
    // because modifying container while iterating it may cause unexpected behavior
//    newEdges.push_back(newE1);
//    newEdges.push_back(newE2);
//    newEdges.push_back(newE3);

    _faces.push_back(newFace1);
    _faces.push_back(newFace2);

    _halfedges.push_back(newHe1);
    _halfedges.push_back(newHe2);
    _halfedges.push_back(newHe3);
    _halfedges.push_back(newHe4);
    _halfedges.push_back(newHe5);
    _halfedges.push_back(newHe6);
}

void MeshDataStructure::flip(const Edge *e)
{
    // halfedges in triangle 1
    Halfedge *he1 = e->halfedge;
    Halfedge *he2 = he1->next;
    Halfedge *he3 = he2->next;

    // halfedges in triangle 2
    Halfedge *he4 = e->halfedge->twin;
    Halfedge *he5 = he4->next;
    Halfedge *he6 = he5->next;

    Vertex *va = he1->vertex;
    Vertex *vb = he2->vertex;
    Vertex *vc = he3->vertex;
    Vertex *vd = he6->vertex;

    Face *f1 = he1->face;
    Face *f2 = he4->face;

    he1->vertex = vd;
    he4->vertex = vc;

    va->halfedge = he5;
    vb->halfedge = he2;

    he1->next = he3;
    he2->next = he4;
    he3->next = he5;
    he4->next = he6;
    he5->next = he1;
    he6->next = he2;

    he1->face = f1;
    he3->face = f1;
    he5->face = f1;

    he2->face = f2;
    he4->face = f2;
    he6->face = f2;

    f1->halfedge = he1;
    f2->halfedge = he4;
}

bool MeshDataStructure::canFlip(const Edge *e)
{
    Vertex *v1 = e->halfedge->vertex;
    Vertex *v2 = e->halfedge->twin->vertex;

    unordered_set<Vertex*> neighbors1 = getNeighborVertices(v1);
    unordered_set<Vertex*> neighbors2 = getNeighborVertices(v2);

    if (neighbors1.size() == 3 || neighbors2.size() == 3) {
        return false;
    }
    return true;
}

int getDegreeDeviation(int degree, int standard)
{
    return abs(degree - standard);
}

bool MeshDataStructure::validToFlipForRemesh(Edge *e)
{
    // halfedges in triangle 1
    Halfedge *he1 = e->halfedge;
    Halfedge *he2 = he1->next;
    Halfedge *he3 = he2->next;

    // halfedges in triangle 2
    Halfedge *he4 = e->halfedge->twin;
    Halfedge *he5 = he4->next;
    Halfedge *he6 = he5->next;

    Vertex *va = he1->vertex;
    Vertex *vb = he2->vertex;
    Vertex *vc = he3->vertex;
    Vertex *vd = he6->vertex;

    int degreeA = getNeighborVertices(va).size();
    int degreeB = getNeighborVertices(vb).size();
    int degreeC = getNeighborVertices(vc).size();
    int degreeD = getNeighborVertices(vd).size();

    int currentTotalDeviation = getDegreeDeviation(degreeA, 6) +
                                getDegreeDeviation(degreeB, 6) +
                                getDegreeDeviation(degreeC, 6) +
                                getDegreeDeviation(degreeD, 6);

    int newDegreeA = degreeA - 1;
    int newDegreeB = degreeB - 1;
    int newDegreeC = degreeC + 1;
    int newDegreeD = degreeD + 1;

    int flippedTotalDeviation = getDegreeDeviation(newDegreeA, 6) +
                                getDegreeDeviation(newDegreeB, 6) +
                                getDegreeDeviation(newDegreeC, 6) +
                                getDegreeDeviation(newDegreeD, 6);

    return flippedTotalDeviation < currentTotalDeviation;
}

void MeshDataStructure::moveVertexTangentially(Vertex *v, Vector3f centroid, float weight)
{
    Vector3f diff = centroid - v->point;
    Vector3f N = getVertexNormal(v);

    // project diff onto tangent plane
    diff = diff - N.dot(diff) * N;

    v->point = v->point + weight * diff;
}

Vector3f MeshDataStructure::getCentroid(Vertex *v)
{
    Vector3f sum(0, 0, 0);
    unordered_set<Vertex*> neighbors = getNeighborVertices(v);

    for (Vertex *q : neighbors) {
        sum += q->point;
    }
    return sum / neighbors.size();
}

Vector3f MeshDataStructure::getVertexNormal(const Vertex *v)
{
    Vector3f sum(0, 0, 0);
    Halfedge *curr = v->halfedge;
    int cnt = 0;

    // traverse neighboring faces of v
    // and sum their normals
    do {
        vector<Vertex*> vertices = getVerticesInTriangle(curr);
        Vector3f faceNormal = getFaceNormal(vertices.at(0), vertices.at(1), vertices.at(2));
        sum += faceNormal;

        curr = curr->twin->next;
        cnt++;
    }
    while (curr != v->halfedge);

    // average the nromals from neighbor faces
    return (sum / cnt).normalized();
}

vector<MeshDataStructure::Vertex*> MeshDataStructure::getVerticesInTriangle(Halfedge *he)
{
    vector<Vertex*> vertices;
    Halfedge *curr = he;

    do {
        vertices.push_back(curr->vertex);
        curr = curr->next;
    }
    while (curr != he);

    return vertices;
}

vector<Vector3f> MeshDataStructure::getVertices()
{
    vector<Vector3f> ret;
    int id = 0;

    for (Vertex *v : _vertices) {
        ret.push_back(Vector3f(v->point));
        _vertexIdMap[v] = id++;
    }
    return ret;
}

// call getVertices() before using this function!
vector<Vector3i> MeshDataStructure::getFaces()
{
    vector<Vector3i> ret;
    for (Face *f : _faces) {
        Vector3i points(0, 0, 0);
        Halfedge *curr = f->halfedge;
        int i = 0;

        do {
            if (!_vertexIdMap.contains(curr->vertex)) {
                cout << "data corrupted, abort" << endl;
                return ret;
            }

            points[i++] = _vertexIdMap[curr->vertex];
            curr = curr->next;
        }
        while (curr != f->halfedge);

        ret.push_back(points);
    }
    return ret;
}
