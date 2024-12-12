CXX = g++
#CXXFLAGS = -std=c++17 -O2 -fopenmp
CXXFLAGS = -std=c++17 -O2 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp
EXEC = multipole
SRC = main.cpp
PYTHON_SCRIPT = plot_trajectories.py
CSV_FILE = particle_data.csv

all: $(EXEC) run_plot

$(EXEC): $(SRC)
	@echo "Compiling C++ program..."
	$(CXX) $(CXXFLAGS) $(SRC) -o $(EXEC)
	@echo "Compilation successful."

run_simulation: $(EXEC)
	@echo "Running C++ simulation program..."
	./$(EXEC)
	@echo "C++ program execution completed."

run_plot: run_simulation
	@echo "Running Python visualization script..."
	python $(PYTHON_SCRIPT)
	@echo "Python script execution completed."

run: all

clean:
	@echo "Cleaning up generated files..."
	rm -f $(EXEC)
	rm -f $(CSV_FILE)
	@echo "Cleanup completed."

.PHONY: all run_simulation run_plot run clean

