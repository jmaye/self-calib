cmake_minimum_required(VERSION 2.6)

include(ReMake)

remake_set(CMAKE_CXX_FLAGS -std=c++0x CACHE FORCE INIT)

remake_project(
  aslam
  VERSION 0.1
  RELEASE alpha
  SUMMARY "A set of data structures and algorithms for robotics"
  AUTHOR "Paul Timothy Furgale"
  AUTHOR "Jerome Maye"
  CONTACT "jerome.maye@mavt.ethz.ch"
  HOME "http://aslforge.ethz.ch/projects/aslam"
  LICENSE "GNU Lesser General Public License (LGPL)"
  INSTALL /usr
)

remake_pack_deb()

remake_pack_deb(
  COMPONENT dev
  DESCRIPTION "Development headers"
  aslam
)

remake_doc(
  html
  latex
)

remake_pack_deb(
  COMPONENT doc
  DESCRIPTION "Documentation"
  aslam
)
