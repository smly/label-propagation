CC   = g++
OBJS = lprop.o graph.o
TARGET = lprop
LIBS = -pthread
INCLUDES = -I.
$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(LIBS) $(INCLUDES) $(OBJS)
.c.o:
	$(CC) -c $<
clean:
	rm $(OBJS) $(TARGET)
