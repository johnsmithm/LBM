CXX = g++
CXXFLAGS = -Wall -pedantic -std=c++11
EXTRA = -ftree-vectorizer-verbose=1 
TARGET1 = lbm
TARGET2 = GrayScaleImage
TARGET3 = lodepng
#HXX= LBM.h  Timer.h imageClass/GrayScaleImage.h

OBJS = $(TARGET1).o ./imageClass/$(TARGET3).o ./imageClass/$(TARGET2).o

EXEC = lbm

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC) 
#all: $(TARGET1)

#$(TARGET1): $(OBJS) $(HXX) 
#	$(CXX) $(CXXFLAGS) -o $(TARGET1) $(OBJS) $(LDFLAGS) $(LIBS)
#%.o: %.cpp
#	$(CXX) $(CXXFLAGS) -c $<

#$(TARGET2): $(OBJS) $(HXX) 
#	$(CXX) $(CXXFLAGS) -o $(TARGET2) $(OBJS) $(LDFLAGS) $(LIBS)
#%.o: %.cpp
#	$(CXX) $(CXXFLAGS) -c $<

$(TARGET1).o: $(TARGET1).cpp
	$(CXX) -c $(CXXFLAGS)  $(TARGET1).cpp  -o $(TARGET1).o


$(TARGET2).o: ./imageClass/$(TARGET3).cpp
	$(CXX) -c $(CXXFLAGS)  ./imageClass/$(TARGET3).cpp  -o $(TARGET3).o

$(TARGET3).o: ./imageClass/$(TARGET2).cpp
	$(CXX) -c $(CXXFLAGS)  ./imageClass/$(TARGET2).cpp  -o $(TARGET2).o

clean:
	rm -rf *.o $(TARGET1)

run1:
	 ./$(EXEC) scenario1	
run2:
	 ./$(EXEC) scenario2