################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Algo_genet.cpp 

OBJS += \
./Algo_genet.o 

CPP_DEPS += \
./Algo_genet.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/courtin/Documents/M2/ProjetC/boost_1_68_0 -I/home/courtin/Documents/M2/ProjetC/detection-of-epistasia-pattern/Genetic/includes -O0 -g3 -Wall -c -fmessage-length=0 -O2 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


