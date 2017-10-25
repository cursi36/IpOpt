cmake_minimum_required(VERSION 2.6)
project(my_solver)

add_executable(my_solver src/main.cpp)

add_executable(main_solver src/main_solver.cpp)

link_directories(/home/francesco/IpOpt/CoinIpopt/Ipopt/src)

include_directories(/home/francesco/IpOpt/CoinIpopt/build/include/coin)

include_directories(include)


install(TARGETS my_solver RUNTIME DESTINATION bin)
