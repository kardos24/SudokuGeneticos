CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -fpermissive

OBJS =		pruebas.o

LIBS =

TARGET =	pruebas.exe

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
