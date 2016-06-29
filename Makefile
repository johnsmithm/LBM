CXX = g++
CXXFLAGS = -Wall -pedantic -std=c++11
EXTRA = -ftree-vectorizer-verbose=1 
TARGET =lbm
HXX= LBM.h  Timer.h

OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(OBJS) $(HXX) 
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf *.o $(TARGET)
