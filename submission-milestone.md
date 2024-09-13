## Mesh (milestone submission)

Please fill this out and submit your work to Gradescope by the milestone deadline.

### Mesh Validator
Describe what your mesh validator checks for here. This can be a list of assertions.

My mesh data structure adopts halfedge to represent geometries.
Each face refers to a halfedge in the face.
Each edge refers to a halfedge in the edge.
Each vertex refers to a halfedge connecting the vertex.
Each halfedge refers to a starting vertex, a twin halfedge, a next halfedge, an edge it belongs to, and a face it belongs to.

My validator validates the following:
- validates each halfedge is not missing data
- validates each edge is not missing data
- validates each face is not missing data
- validates each vertex is not missing data
- validates each halfedge's twin's twin points back to itself
- validates each halfedge's vertex and that of its twin is different
- validates each vertex's halfedge's vertex points back to itself
- traverses each halfedge's next to validate it points back to itself
- detects if there is loop when traversing each halfedge's next
- traverses halfedges around each vertex to validate it points back to itself
- detects if there is loop when traversing halfedges around each vertex

### Collaboration/References

### Known Bugs
