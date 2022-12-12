

using namespace globalVariable;

TEST_CASE("Test class Mesh 2D", "[mesh]") {

   Mesh2 mesh(3, 3, 0, 0, 2, 2);

   REQUIRE(mesh.getNbElement() == 8);
   REQUIRE(mesh.getNbNode() == 9);
   REQUIRE(mesh.getNbBorder() == 8);

   const auto &T(mesh[0]); // unit triangle
   REQUIRE(mesh(T) == 0);
   REQUIRE(isEqual(T[0][0], 0.));
   REQUIRE(isEqual(T[0][1], 0.));
   REQUIRE(isEqual(T[1][0], 1.));
   REQUIRE(isEqual(T[1][1], 0.));
   REQUIRE(isEqual(T[2][0], 0.));
   REQUIRE(isEqual(T[2][1], 1.));

   const auto &V(mesh(4));
   REQUIRE(mesh(V) == 4);
   REQUIRE(isEqual(V[0], 1.));
   REQUIRE(isEqual(V[1], 1.));

   const auto &BE(mesh.be(1));
   REQUIRE(mesh(BE) == 1);
   REQUIRE(isEqual(BE[0][0], 1.));
   REQUIRE(isEqual(BE[0][1], 0.));
   REQUIRE(isEqual(BE[1][0], 2.));
   REQUIRE(isEqual(BE[1][1], 0.));
   REQUIRE(BE.lab);

   // Check the index of the edge of the triangle 0
   typedef SortArray<int, 2> N2;
   REQUIRE(mesh.itemadj(0, 0) == N2(1, 3));
   REQUIRE(mesh.itemadj(0, 1) == N2(0, 3));
   REQUIRE(mesh.itemadj(0, 2) == N2(0, 1));

   // Check the neighbor element (first element in the tuple. The second element
   // is the face idx in the neighbor element)
   REQUIRE(std::get<0>(mesh.getElementAdj(0, 1)) == -1); // -1 is border element
   REQUIRE(std::get<0>(mesh.getElementAdj(0, 2)) == -1);
   REQUIRE(std::get<0>(mesh.getElementAdj(3, 1)) == -1);

   REQUIRE(mesh.getElementAdj(0, 0) == std::tuple<int, int>{1, 0});
   REQUIRE(mesh.getElementAdj(3, 0) == std::tuple<int, int>{2, 0});

   // get the element and the local idx of the face corresponding to the
   // boundary element
   REQUIRE(mesh.getBoundaryElement(0) == std::tuple<int, int>{0, 2});
   REQUIRE(mesh.getBoundaryElement(1) == std::tuple<int, int>{2, 2});
   REQUIRE(mesh.getBoundaryElement(2) == std::tuple<int, int>{3, 1});
   REQUIRE(mesh.getBoundaryElement(3) == std::tuple<int, int>{7, 1});
   REQUIRE(mesh.getBoundaryElement(4) == std::tuple<int, int>{5, 2});
   REQUIRE(mesh.getBoundaryElement(5) == std::tuple<int, int>{7, 2});
   REQUIRE(mesh.getBoundaryElement(6) == std::tuple<int, int>{0, 1});
   REQUIRE(mesh.getBoundaryElement(7) == std::tuple<int, int>{4, 1});
}

TEST_CASE("Test class MeshQuad 2D", "[mesh]") {

   MeshQuad2 mesh(3, 3, 0, 0, 2, 2);

   REQUIRE(mesh.getNbElement() == 4);
   REQUIRE(mesh.getNbNode() == 9);
   REQUIRE(mesh.getNbBorder() == 8);

   const auto &T(mesh[0]); // first element and the unit square
   REQUIRE(mesh(T) == 0);  // index of element T
   // T[i] returns R2 coordinates of vertex with index i of triangle T
   REQUIRE(isEqual(T[0][0], 0.));
   REQUIRE(isEqual(T[0][1], 0.));
   REQUIRE(isEqual(T[1][0], 1.));
   REQUIRE(isEqual(T[1][1], 0.));
   REQUIRE(isEqual(T[2][0], 1.));
   REQUIRE(isEqual(T[2][1], 1.));
   REQUIRE(isEqual(T[3][0], 0.));
   REQUIRE(isEqual(T[3][1], 1.));

   const auto &V(mesh(4)); // vertex with index 4
   REQUIRE(mesh(V) == 4);
   REQUIRE(isEqual(V[0], 1.));
   REQUIRE(isEqual(V[1], 1.));

   const auto &BE(mesh.be(2)); // border element with index 2
   REQUIRE(mesh(BE) == 2);
   // BE[i] returns coordinates of vertex i of the border element
   REQUIRE(isEqual(BE[0][0], 2.));
   REQUIRE(isEqual(BE[0][1], 0.));
   REQUIRE(isEqual(BE[1][0], 2.));
   REQUIRE(isEqual(BE[1][1], 1.));
   REQUIRE(BE.lab == 2);

   // Check the vertex indices of the edges of elements
   typedef SortArray<int, 2> N2;
   // Input: first index = element, second index = local face index
   // first square
   REQUIRE(mesh.itemadj(0, 0) == N2(0, 1));
   REQUIRE(mesh.itemadj(0, 1) == N2(1, 4));
   REQUIRE(mesh.itemadj(0, 2) == N2(3, 4));
   REQUIRE(mesh.itemadj(0, 3) == N2(0, 3));
   // last square
   REQUIRE(mesh.itemadj(3, 0) == N2(4, 5));
   REQUIRE(mesh.itemadj(3, 1) == N2(5, 8));
   REQUIRE(mesh.itemadj(3, 2) == N2(7, 8));
   REQUIRE(mesh.itemadj(3, 3) == N2(4, 7));

   // Check the neighbor element (first element in the tuple. The second element
   // is the face idx in the neighbor element)
   REQUIRE(std::get<0>(mesh.getElementAdj(0, 0)) == -1); // -1 is border element
   REQUIRE(std::get<0>(mesh.getElementAdj(0, 3)) == -1);
   REQUIRE(std::get<0>(mesh.getElementAdj(3, 1)) == -1);
   REQUIRE(std::get<0>(mesh.getElementAdj(3, 2)) == -1);

   // Check that element and local face indexing agrees between adjacent
   // elements
   REQUIRE(mesh.getElementAdj(0, 1) == std::tuple<int, int>{1, 3});
   REQUIRE(mesh.getElementAdj(1, 3) == std::tuple<int, int>{0, 1});
   REQUIRE(mesh.getElementAdj(0, 2) == std::tuple<int, int>{2, 0});
   REQUIRE(mesh.getElementAdj(2, 0) == std::tuple<int, int>{0, 2});
   REQUIRE(mesh.getElementAdj(2, 1) == std::tuple<int, int>{3, 3});
   REQUIRE(mesh.getElementAdj(3, 3) == std::tuple<int, int>{2, 1});

   // Get the element and the local idx of the face corresponding to the
   // boundary element
   REQUIRE(mesh.getBoundaryElement(0) == std::tuple<int, int>{0, 0});
   REQUIRE(mesh.getBoundaryElement(1) == std::tuple<int, int>{1, 0});
   REQUIRE(mesh.getBoundaryElement(2) == std::tuple<int, int>{1, 1});
   REQUIRE(mesh.getBoundaryElement(3) == std::tuple<int, int>{3, 1});
   REQUIRE(mesh.getBoundaryElement(4) == std::tuple<int, int>{2, 2});
   REQUIRE(mesh.getBoundaryElement(5) == std::tuple<int, int>{3, 2});
   REQUIRE(mesh.getBoundaryElement(6) == std::tuple<int, int>{0, 3});
   REQUIRE(mesh.getBoundaryElement(7) == std::tuple<int, int>{2, 3});
}
