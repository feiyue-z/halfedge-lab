## Mesh (final submission)

Please fill this out and submit your work to Gradescope by the deadline.

### Output Comparison
Run the program with the specified `.ini` config file to compare your output against the reference images. The program should automatically save the output mesh to the `student_outputs/final` folder. Please take a screenshot of the output mesh and place the image in the table below. Do so by placing the screenshot `.png` in the `student_outputs/final` folder and inserting the path in the table.

- For instance, after running the program with the `subdivide_icosahedron_4.ini` config file, go to and open `student_outputs/final/subdivide_icosahedron_4.obj`. Take a screenshot of the mesh and place the screenshot in the first row of the first table in the column titled `Your Output`.
- The markdown for the row should look something like `| subdivide_icosahedron_4.ini |  ![](ground_truth_pngs/final/subdivide_icosahedron_4.png) | ![](student_outputs/final/subdivide_icosahedron_4.png) |`

If you are not using the Qt framework, you may also produce your outputs otherwise so long as the output images show up in the table. In this case, please also describe how your code can be run to reproduce your outputs.

> Qt Creator users: If your program can't find certain files or you aren't seeing your output images appear, make sure to:<br/>
> 1. Set your working directory to the project directory
> 2. Set the command-line argument in Qt Creator to `template_inis/final/<ini_file_name>.ini`

Note that your outputs do **not** need to exactly match the reference outputs. There are several factors that may result in minor differences, especially for certain methods like simplification where equal-cost edges may be handled differently.



Please do not attempt to duplicate the given reference images; we have tools to detect this.

| `.ini` File To Produce Output | Expected Output | Your Output |
| :---------------------------------------: | :--------------------------------------------------: | :-------------------------------------------------: | 
| subdivide_icosahedron_4.ini |  ![](ground_truth_pngs/final/subdivide_icosahedron_4.png) | ![Place screenshot of student_outputs/final/subdivide_icosahedron_4.obj here](student_outputs/final/subdivide_icosahedron_4.png) |
| simplify_sphere_full.ini |  ![](ground_truth_pngs/final/simplify_sphere_full.png) | ![Place screenshot of student_outputs/final/simplify_sphere_full.obj here](student_outputs/final/simplify_sphere_full.png) |
| simplify_cow.ini | ![](ground_truth_pngs/final/simplify_cow.png) | ![Place screenshot of student_outputs/final/simplify_cow.obj here](student_outputs/final/simplify_cow.png) |

Output for Isotropic Remeshing (Note: if you did not implement this you can just skip this part)
| `.ini` File To Produce Output | Input Mesh .png | Remeshed Mesh .png |
| :---------------------------------------: | :--------------------------------------------------: | :-------------------------------------------------: | 
| <template_inis/final/remesh_peter.ini> |  ![Place screenshot of your input mesh here](student_outputs/final/peter.png) | ![Place screenshot of your remeshed mesh here](student_outputs/final/remeshed_peter.png) |



Output for Bilateral Mesh Denoising (Note: if you did not implement this you can just skip this part)
| `.ini` File To Produce Output | Noisy Mesh .png | Denoised Mesh .png |
| :---------------------------------------: | :--------------------------------------------------: | :-------------------------------------------------: | 
| <template_inis/final/denoise_bunny.ini> |  ![Place screenshot of a noisy mesh here](student_outputs/final/noised_bunny.png) | ![Place screenshot of your denoised mesh here](student_outputs/final/denoised_bunny.png) |



Output for any other Geometry Processing Functions (Note: if you did not implement this you can just skip this part)
| `.ini` File To Produce Output | Input | Output |
| :---------------------------------------: | :--------------------------------------------------: | :-------------------------------------------------: | 
| <Path to your .ini file> |  ![Place screenshot input mesh here]() | ![Place screenshot of output mesh here]() |


### Design Choices

#### Mesh Data Structure 
Describe your mesh data structure here. 

My mesh data structure adopts halfedge to represent geometries.
Each face refers to a halfedge in the face.
Each edge refers to a halfedge in the edge.
Each vertex refers to a halfedge connecting the vertex.
Each halfedge refers to a starting vertex, a twin halfedge, a next halfedge, an edge it belongs to, and a face it belongs to.

#### Mesh Validator
Describe what your mesh validator checks for here. This can be a list.

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

#### Run Time/Efficency 
Describe how you achieved efficient asymptotic running times for your geometry processing functions, including the data structures you used.

The data structures I used include Map, Vector, Multiset, and Unordered Set. I run through mesh data once without nested traversal, so the time complexity should be linear at maximum.

### Collaboration/References

### Known Bugs
