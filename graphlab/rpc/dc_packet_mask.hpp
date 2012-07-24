/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#ifndef DC_PACKET_MASK_HPP
#define DC_PACKET_MASK_HPP
namespace graphlab {
  // ---------  Packet header types --------------

  /**
   * \internal
   * \ingroup rpc
   * Used for regular calls which go into a thread pool
   * for evaluation
   */
  const unsigned char STANDARD_CALL = 1;


  /**
   * \internal
   * \ingroup rpc
   * 
   * If WAIT_FOR_REPLY is set, the function call's
  return will be passed back to the caller */
  const unsigned char WAIT_FOR_REPLY = 4;

  /**
   * \internal
    \ingroup rpc
   * 
    If control packet flag is set, this packet 
    does not increment any counters.
  */
  const unsigned char CONTROL_PACKET = 16; 
  
  /**
   * \internal
   * \ingroup rpc
   * 
   * Used to identify that this packet was 
   * a reply to a previous request.
   */
  const unsigned char REPLY_PACKET = 32;

  /**
   * \internal
   * \ingroup rpc
   * 
   * Used to identify that this packet is
   * serialized using a POD mechanism;
   */
  const unsigned char POD_CALL = 64;
}
#endif

