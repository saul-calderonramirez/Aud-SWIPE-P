################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../AuditiveSystem.c \
../Clocks.c \
../MPIMatlabCommunicator.c \
../Matrix.c \
../PrimeMultiF0.c \
../SPUtilities.c \
../TestFactory.c \
../UnitTestsMPI.c \
../main.c 

OBJS += \
./AuditiveSystem.o \
./Clocks.o \
./MPIMatlabCommunicator.o \
./Matrix.o \
./PrimeMultiF0.o \
./SPUtilities.o \
./TestFactory.o \
./UnitTestsMPI.o \
./main.o 

C_DEPS += \
./AuditiveSystem.d \
./Clocks.d \
./MPIMatlabCommunicator.d \
./Matrix.d \
./PrimeMultiF0.d \
./SPUtilities.d \
./TestFactory.d \
./UnitTestsMPI.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


