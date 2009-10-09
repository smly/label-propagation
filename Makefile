CC   = g++
OBJS = main.o lprop.o graph.o
TARGET = lprop
LIBS = -pthread
INCLUDES = -I.
$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(LIBS) $(INCLUDES) $(OBJS)
.c.o:
	$(CC) -c $<
clean:
	rm $(OBJS) $(TARGET)
