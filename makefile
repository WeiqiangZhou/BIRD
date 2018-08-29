all:
	g++ -std=c++11 source/predict_main.cpp source/predict_function.cpp -o BIRD_predict
	g++ source/create_library_main.cpp source/create_library_function.cpp -o BIRD_build_library
