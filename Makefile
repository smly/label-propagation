CC   = g++
#OBJS = dwalk.o utils.o graph.o
#TARGET = dwalk
OBJS = label_propagation.o utils.o graph.o
TARGET = label_propagation
LIBS = -lboost_filesystem -lglog -pthread -lgflags
INCLUDES = -I.
$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(LIBS) $(INCLUDES) $(OBJS)
.c.o:
	$(CC) -c $<
clean:
	rm $(OBJS) $(TARGET)