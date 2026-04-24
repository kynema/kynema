#pragma once

namespace kynema_fmb::util {

extern "C" {

/// Implement a test controller that returns pitch angle in radians (ranges from -90 to 90 starting
/// at zero) - used for testing purposes
void PITCH_CONTROLLER(
    float avrSWAP[], int* aviFAIL, const char* accINFILE, char* avcOUTNAME, char* avcMSG
);

}  // extern "C"

}  // namespace kynema_fmb::util
