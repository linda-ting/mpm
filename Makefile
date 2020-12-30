CPP = g++ -std=c++11
CPPFLAGS = -Wall -DNDEBUG # -g
LDFLAGS = -lm

all: material_point

material_point: main.cpp libmesh_query.a
	$(CPP) $(CPPFLAGS) -I /usr/local/include/eigen3 main.cpp libmesh_query.a -o material_point

libmesh_query.a: mesh_query.o predicates.o bounding_box_tree.o
	ar r $@ mesh_query.o predicates.o bounding_box_tree.o;
	ranlib $@

mesh_query.o: mesh_query.cpp mesh_query.h bounding_box_tree.h bounding_box.h predicates.h vec.h util.h
	$(CPP) $(CPPFLAGS) -o $@ -c mesh_query.cpp

predicates.o: predicates.cpp predicates.h
	$(CPP) $(CPPFLAGS) -o $@ -c predicates.cpp

bounding_box_tree.o: bounding_box_tree.cpp bounding_box_tree.h bounding_box.h vec.h util.h
	$(CPP) $(CPPFLAGS) -o $@ -c bounding_box_tree.cpp

clean:
	-rm material_point libmesh_query.a mesh_query.o predicates.o bounding_box_tree.o
