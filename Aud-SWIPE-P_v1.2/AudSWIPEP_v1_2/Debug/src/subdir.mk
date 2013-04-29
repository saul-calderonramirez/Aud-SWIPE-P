################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/AudSWIPEP.c \
../src/AuditiveSystem.c \
../src/Clocks.c \
../src/MPIMatlabCommunicator.c \
../src/Matrix.c \
../src/SPUtilities.c \
../src/TestFactory.c \
../src/UnitTestsMPI.c \
../src/main.c 

OBJS += \
./src/AudSWIPEP.o \
./src/AuditiveSystem.o \
./src/Clocks.o \
./src/MPIMatlabCommunicator.o \
./src/Matrix.o \
./src/SPUtilities.o \
./src/TestFactory.o \
./src/UnitTestsMPI.o \
./src/main.o 

C_DEPS += \
./src/AudSWIPEP.d \
./src/AuditiveSystem.d \
./src/Clocks.d \
./src/MPIMatlabCommunicator.d \
./src/Matrix.d \
./src/SPUtilities.d \
./src/TestFactory.d \
./src/UnitTestsMPI.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I/usr/include/mpi -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


