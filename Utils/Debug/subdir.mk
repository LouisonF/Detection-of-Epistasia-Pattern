################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Miscellaneous.cpp \
../Parametersfileparsing.cpp \
../datainput.cpp \
../main_test.cpp 

OBJS += \
./Miscellaneous.o \
./Parametersfileparsing.o \
./datainput.o \
./main_test.o 

CPP_DEPS += \
./Miscellaneous.d \
./Parametersfileparsing.d \
./datainput.d \
./main_test.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/Utils/includes -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


