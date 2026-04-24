#pragma once

namespace kynema_fmb::util {

extern "C" {

/// @brief This Bladed-style DLL controller is used to implement a variable-speed generator-torque
/// controller and PI collective blade pitch controller for the NREL Offshore 5MW baseline wind
/// turbine. This routine was originally written in Fortran by J. Jonkman of NREL/NWTC for use in the
/// IEA Annex XXIII OC3 studies.
/// @param avrSWAP The swap array, used to pass data to, and receive data from, the DLL
///                controller.
/// @param aviFAIL Flag used to indicate the success of this DLL call set as follows:
///                0 if the DLL call was successful,
///                > 0 if the DLL call was successful but cMessage should be issued as a warning
///                messsage,
///                < 0 if the DLL call was unsuccessful or for any other reason the simulation
///                is to be stopped at this point with cMessage as the error message
/// @param accINFILE The name of the parameter input file, 'DISCON.IN'
/// @param avcOUTNAME OUTNAME (Simulation RootName)
/// @param avcMSG MESSAGE (Message from DLL to simulation code [ErrMsg])  The message which will
///               be displayed by the calling program if aviFAIL <> 0
void DISCON(float avrSWAP[], int* aviFAIL, const char* accINFILE, char* avcOUTNAME, char* avcMSG);

}  // extern "C"

}  // namespace kynema_fmb::util
