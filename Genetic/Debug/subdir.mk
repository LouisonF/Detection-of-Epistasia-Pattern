################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Algo_genet.cpp \
../Child.cpp \
../ContingencyTable.cpp \
../G2test.cpp \
../InitialMatrix.cpp \
../Parent.cpp \
../Population.cpp \
../TheoricalTable.cpp \
../old_algo.cpp \
../test.cpp 

OBJS += \
./Algo_genet.o \
./Child.o \
./ContingencyTable.o \
./G2test.o \
./InitialMatrix.o \
./Parent.o \
./Population.o \
./TheoricalTable.o \
./old_algo.o \
./test.o 

CPP_DEPS += \
./Algo_genet.d \
./Child.d \
./ContingencyTable.d \
./G2test.d \
./InitialMatrix.d \
./Parent.d \
./Population.d \
./TheoricalTable.d \
./old_algo.d \
./test.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/courtin/Documents/M2/ProjetC/boost_1_68_0 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

