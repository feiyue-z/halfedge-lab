#pragma once

#include <map>
#include <queue>
#include <set>
#include <unordered_set>

#include "mesh.h"

typedef std::pair<Eigen::Vector3f, Eigen::Vector3f> PointPair;

class MeshDataStructure
{
public:
    struct Vertex;
    struct Edge;
    struct Face;

    struct Halfedge
    {
        Halfedge *twin = nullptr;
        Halfedge *next = nullptr;
        Vertex *vertex = nullptr;
        Edge *edge = nullptr;
        Face *face = nullptr;
    };

    struct Vertex
    {
        Halfedge *halfedge = nullptr;
        Eigen::Vector3f point = Eigen::Vector3f(0, 0, 0);

        Vertex(Eigen::Vector3f p) : point(p) {}
    };

    struct Edge
    {
        Halfedge *halfedge;
        float cost = 0.f;
        Vertex *optimalVertex;
    };

    struct Face
    {
        Halfedge *halfedge;
    };

    struct Vector3fComparator {
        bool operator()(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2) const {
            if (v1.x() < v2.x()) return true;
            if (v1.x() > v2.x()) return false;
            if (v1.y() < v2.y()) return true;
            if (v1.y() > v2.y()) return false;
            return v1.z() < v2.z();
        }
    };

    struct PointPairComparator {
        bool operator()(const PointPair& o1, const PointPair& o2) const {
            Vector3fComparator comp;

            auto sorted_o1 = std::minmax(o1.first, o1.second, comp);
            auto sorted_o2 = std::minmax(o2.first, o2.second, comp);

            if (comp(sorted_o1.first, sorted_o2.first)) return true;
            if (comp(sorted_o2.first, sorted_o1.first)) return false;

            return comp(sorted_o1.second, sorted_o2.second);
        }
    };

    struct EdgeCostComparator {
        bool operator()(const Edge *o1, const Edge *o2) const {
            return o1->cost < o2->cost;
        }
    };

    void init(const Mesh &mesh);
    void validate();
    void subdivide(int numIterations);
    void simplify(int numFaceRemoves);
    void noise(int numIterations, float scale);
    void denoise(int numIterations, float sigma_c, float sigma_s, float rho);
    void remesh(int numIterations, float weight);

    std::vector<Eigen::Vector3f> getVertices();
    std::vector<Eigen::Vector3i> getFaces();

private:
    std::vector<Vertex*> _vertices;
    std::vector<Halfedge*> _halfedges;
    std::vector<Edge*> _edges;
    std::vector<Face*> _faces;
    std::map<Vertex*, int> _vertexIdMap;

    // atomic
    void splitEdge(Edge *e,
                   std::map<Edge*, Vertex*> &newVertexMap);
    void flip(const Edge *e);
    void collapse(Edge *eCollap);
    void collapseAndUpdateCost(Edge *e,
                               std::map<Vertex*, Eigen::Matrix4f> &vertexQuadricMap,
                               std::multiset<Edge*, EdgeCostComparator> &edgeCostSet);

    void subdivideOnce();
    void updateStructure(const std::vector<Vertex*> &faceVertices,
                         std::vector<Halfedge*> &halfedges,
                         std::vector<Edge*> &edges,
                         std::vector<Face*> &faces,
                         std::map<PointPair, Halfedge*, PointPairComparator> &halfedgeMap);
    void computeVertexQuadric(const Face *f,
                              std::map<Vertex*, Eigen::Matrix4f> &vertexQuadricMap);
    void mergeCollapsedEdges(const Halfedge *he);
    void ensureValidHalfedgeForVertex(Vertex *v,
                                      std::unordered_set<Halfedge*> heRemove);
    void eraseEdgeFromSet(Edge *e,
                          std::multiset<Edge*, EdgeCostComparator> &edgeCostSet);
    void moveVertexTangentially(Vertex *v,
                                Eigen::Vector3f centroid,
                                float weight);
    void splitTwoTrianglesIntoFour(Edge *e,
                                   Vertex *vMid);

    // checker
    bool canFlip(const Edge *e);
    bool canCollapse(Edge *e);
    bool edgeIsInSet(Edge *e,
                     std::multiset<Edge*, EdgeCostComparator> &edgeCostSet);
    bool validToFlipForRemesh(Edge *e);

    float getEdgeCostAndOptimalVertex(Edge *e,
                                      std::map<Vertex*, Eigen::Matrix4f> &vertexQuadricMap);
    float getMeanEdgeLength();
    float getEdgeLength(Edge *e);

    Eigen::Vector3f getVertexNormal(const Vertex *v);
    Eigen::Vector3f getFaceNormal(Vertex *v1, Vertex *v2, Vertex *v3);
    Eigen::Vector3f getDenoisedPoint(const Vertex* v,
                                    float sigma_c,
                                    float sigma_s,
                                    float rho);
    Eigen::Vector3f getCentroid(Vertex *v);

    std::vector<Vertex*> getVerticesInTriangle(Halfedge *he);
    std::unordered_set<Vertex*> getNeighborVertices(const Vertex *v);
    std::unordered_set<Halfedge*> getHalfedgesInTriangle(Halfedge *he);
    std::unordered_set<Edge*> updateCollapsedVertex(const Edge* e,
                                                    const Vertex* vOld,
                                                    Vertex* vNew);
};
