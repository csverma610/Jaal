OBJS = testPoisson2D.o

OBJS += ../AllQuadMeshGenerator.o
OBJS += ../Mesh.o
OBJS += ../MeshTopology.o
OBJS += ../MeshGeometry.o
OBJS += ../MinDiam.o
OBJS += ../MinBox.o
OBJS += ../MeshImporter.o
OBJS += ../MeshExporter.o
OBJS += ../MeshQuality.o
OBJS += ../Logger.o
OBJS += ../MeshRelManager.o
OBJS += ../off_io.o
OBJS += ../obj_io.o
OBJS += ../xml_io.o
OBJS += ../vtk_io.o
OBJS += ../tri_io.o
OBJS += ../gmsh_io.o
OBJS += ../povray_io.o
OBJS += ../medit_io.o
OBJS += ../Node.o
OBJS += ../Edge.o
OBJS += ../Face.o
OBJS += ../Cell.o
OBJS += ../Simplex.o
OBJS += ../Geometry.o
OBJS += ../MeshRefine3D.o
OBJS += ../tfiblend.o
OBJS += ../GeomPredicates.o
OBJS += ../DelaunayMesh2D.o
OBJS += ../BoundingBox.o
OBJS += ../Geometry.o

CPPFLAGS =  -O3 -fopenmp -std=c++17
CPPFLAGS += -I../../include
CPPFLAGS += -I$(EIGEN_DIR)
CPPFLAGS += -I$(MESQUITE_DIR)/include
CPPFLAGS += -I$(TRIANGLE_DIR)/include
CPPFLAGS += -I$(VERDICT_DIR)/include
CPPFLAGS += -I$(ANN_DIR)/include

LIBS     = -lgomp -pthread
LIBS    += -L$(VERDICT_DIR)/lib -lverdict
LIBS    += -L$(ANN_DIR)/lib -lANN

fire:$(OBJS)
	g++ -o fire $(OBJS) $(LIBS)

.o:.cpp
	g++ $(CPPFLAGS) $<

clean:
	\rm -rf *.o
