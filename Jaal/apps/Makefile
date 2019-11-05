OBJS = testPoisson2D.o

OBJS += ../src/QuadMesh/AllQuadMeshGenerator.o
OBJS += ../src/MeshCore/Mesh.o
OBJS += ../src/MeshCore/MeshTopology.o
OBJS += ../src/MeshCore/MeshGeometry.o
OBJS += ../src/MeshCore/MinDiam.o
OBJS += ../src/MeshCore/MinBox.o
OBJS += ../src/MeshCore/MeshImporter.o
OBJS += ../src/MeshCore/MeshExporter.o
OBJS += ../src/MeshCore/MeshQuality.o
OBJS += ../src/MeshCore/Logger.o
OBJS += ../src/MeshCore/MeshRelManager.o
OBJS += ../src/MeshCore/off_io.o
OBJS += ../src/MeshCore/obj_io.o
OBJS += ../src/MeshCore/xml_io.o
OBJS += ../src/MeshCore/vtk_io.o
OBJS += ../src/MeshCore/tri_io.o
OBJS += ../src/MeshCore/gmsh_io.o
OBJS += ../src/MeshCore/povray_io.o
OBJS += ../src/MeshCore/medit_io.o
OBJS += ../src/MeshCore/Node.o
OBJS += ../src/MeshCore/Edge.o
OBJS += ../src/MeshCore/Face.o
OBJS += ../src/MeshCore/Cell.o
OBJS += ../src/MeshCore/Simplex.o
OBJS += ../src/MeshCore/Geometry.o
OBJS += ../src/MeshCore/MeshRefine3D.o
OBJS += ../src/MeshCore/tfiblend.o
OBJS += ../src/MeshCore/GeomPredicates.o
OBJS += ../src/TriMesh/DelaunayMesh2D.o
OBJS += ../src/MeshCore/BoundingBox.o

CPPFLAGS =  -O3 -fopenmp -std=c++17
CPPFLAGS += -I../include
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
