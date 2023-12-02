# Server side manual

need to `import di_zkp_interface`

The functionalities of a server is encapsulated in the python class `di_zkp_interface.ServerInterface`

A test run is in [test.py](../python/test.py).

## Creating a python class instance
To construct an instance of the server interface, run `di_zkp_interface.ServerInterface(num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
weight_bits, random_normal_bit_shifter,
num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
check_type,
b_precomp)`, where
- `num_clients`: number of clients
- `max_malicious_clients`: maximum number of malicious clients
- `dim`: number of parameters
- `num_blinds_per_group_element`: number of blinds per group element e, see [Precomputation & multiple blinds per group element](precomp.md)
- `weight_bits`: number of bits of weight updates
- `random_normal_bit_shifter`: random normal samples are multiplied by 2^`random_normal_bit_shifter` and rounded to the nearest integer. The paper uses `24`. 
- `num_norm_bound_samples`: number of multidimensional normal samples used in proving the l2 norm bound
- `inner_prod_bound_bits`: the number of bits of each inner product between the model update and discretized multidimensional normal sample. Recommended: use `weight_bits + random_normal_bit_shifter + 4`.
- `max_bound_sq_bits`: the maximum number of bits of the sum of squares of inner products. Recommended: use `2 * (weight_bits + random_normal_bit_shifter) + 20`.
- `check_type`: TODO
- `b_precomp`: boolean value, see [Precomputation & multiple blinds per group element](precomp.md).
  - `True`: store precomputations of powers of group elements for faster commitment
  - `False`: don't store precomputations

[//]: # (- `protocol_type`: the type of protocol invoked)

[//]: # (  - `di_zkp_interface.PROTOCOL_TYPE_PRIV`: the protocol in the paper)

[//]: # (  - `di_zkp_interface.PROTOCOL_TYPE_NON_PRIV_INT`: clients send weight updates in clear text &#40;no ZKP involved&#41;. The probabilistic checking is the same. This is used in experiments that compare accuracies between probabilistic checking and strict checking. )

[//]: # (  - `di_zkp_interface.PROTOCOL_TYPE_NON_PRIV_FLOAT`:)

[//]: # (    clients send weight updates in clear text &#40;no ZKP involved&#41;, AND weight updates are not converted into fixed-point integers, AND in the probabilistic checking, random normal samples are not discretized. )
## Initialize FL training
Given a string `seed`, the class method
`initialize_from_seed(seed)` initializes an FL training process. The string `seed` is used to generate independent group elements w<sub>j</sub>. It should be agreed among the server and all the clients. After initialization,we can start the iterative training process. For every iteration, we should call the following 14 methods:

## Class methods to call within every iteration

### 1. Initializing a new iteration
`initialize_new_iteration(norm_bound)` (no return value) sets the l2 norm bound of every client's weight updates at this round of iteration to `norm_bound`. It should be called exactly once at the start of the iteration. 

### 2. Round 1: Process messages received from clients
`receive_1(bytes_str, i)` (no return value) processes the string `bytes_str` which should be sent by client `i` in round 1. The string `bytes_str` should encode:
- the DHKE public key of client `i`, signed with client `i`'s private key corresponding to its public key on the public bulletin board. (See [Public Bulletin Board](bulletin.md).)

`receive_1(bytes_str, i)` should be called exactly once for every `1 <= i <= num_clients`. (Note: client IDs start from `1`. There is no client `0`.)

### 3. Round 2: Generate a message that should be broadcast to all the clients
`send_2()` returns a string that should be broadcast to all the clients in round 2. The return value encodes:
- the signed DHKE public keys from all the clients. 

It should be called exactly once.

### 4. Round 2: Process messages received from clients
`receive_2(bytes_str, i)` (no return value) processes the string `bytes_str` which should be sent by client `i` in round 2. The string `bytes_str` should encode:
- the commitments of the weight updates of client `i`,
- the check string of client `i`'s blind r<sub>i</sub> and the DHKE-encrypted shares of r<sub>i</sub> via verifiable Shamir sharing. 

`receive_2(bytes_str, i)` should be called exactly once for every `1 <= i <= num_clients`.

### 5. Concurrent process before round 3
`concurrent_process_before_send_3()` (no return value) should be called before round 3. It can be run concurrently as soon as the current round of iteration starts.
It samples the multidimensional normal variables, multiplies them by 2^`random_normal_bit_shifter`, rounds them to the nearest integers and computes the corresponding multi-exponentiations of group elements w<sub>j</sub>. It should be called exactly once.

### 6. Round 3: Generate messages that should be sent to clients
`send_3(i)` returns a string that should be sent to client `i` in round 3. It should be called exactly once for every `1 <= i <= num_clients`. The return value encodes:
- the check string of client `j`'s blind r<sub>j</sub> and  the DHKE-encrypted shares of client `j`'s blind r<sub>j</sub> for every `j != i`, 
- the seed that generates the random normal samples, 
- the corresponding multi-exponentiation values.

`send_3(i)` should be called exactly once for every `1 <= i <= num_clients`.

### 7. Round 3: Process messages received from clients
`receive_3(bytes_str, i)` (no return value) processes the string `bytes_str` which should be sent by client `i` in round 3. The string `bytes_str` should encode:
- the proof that client `i`'s committed weight update has l2 norm at most `norm_bound` by probabilistic checking, 
- the list of clients whose shares do not pass their check strings.

`receive_3(bytes_str, i)` should be called exactly once for every `1 <= i <= num_clients`.

### 8. Process before round 4
`process_before_send_4()` (no return value) should be called after round 3 and before round 4. It checks the all proofs, flags the ones with invalid proofs as malicious, and processes the boolean matrix of whether client `i` reports that client `j`'s share passes the check string or not. It should be called exactly once.

### 9. Round 4: Generate messages that should be sent to clients
`send_4(i)` returns a string that should be sent to client `i` in round 4. It should be called exactly once for every `1 <= i <= num_clients`. The return value encodes:
- the list of `j`s that report that client `i`'s share do not pass the check, provided that the list has length at most `max_malicious_clients`.

`send_4(i)` should be called exactly once for every `1 <= i <= num_clients`.

### 10. Round 4: Process messages received from clients
`receive_4(bytes_str, i)` (no return value) processes the string `bytes_str` which should be sent by client `i` in round 4. The string `bytes_str` should encode:
- the clear text shares of r<sub>i</sub> of the clients in the return value of `send_4(i)`, provided that the list has length at most `max_malicious_clients`.

`receive_4(bytes_str, i)` should be called exactly once for every `1 <= i <= num_clients`.

### 11. Process before round 5
`process_before_send_5()` (no return value) should be called after round 4 and before round 5. It checks all the clear text shares with their check strings and flags more clients based on it. It should be called exactly once. 

### 12. Round 5: Generate messages that should be sent to clients
`send_5(i)` returns a string that should be sent to client `i` in round 5. It should be called exactly once for every `1 <= i <= num_clients`. The return value encodes:
- the clear text values of client `i`'s share of r<sub>j</sub> for `j`'s that don't pass the check strings,
- the list of malicious clients.

`send_5(i)` should be called exactly once for every `1 <= i <= num_clients`.

### 13. Round 5: Process messages received from clients
`receive_5(bytes_str, i)` (no return value) processes the string `bytes_str` which should be sent by client `i` in round 5. The string `bytes_str` should encode:
- the sum of client `i`'s shares of r<sub>j</sub> over all the valid client `j`s.

`receive_5(bytes_str, i)` should be called exactly once for every `1 <= i <= num_clients`.

### 14. Finish iteration
`finish_iteration()` computes the aggregated weight updates of valid clients. The following class members (data type: vector of floats) are updated:
- `final_update_float`: the aggregated weight updates of valid clients
- `final_update_float_avg`: the average weight updates of valid clients

`finish_iteration()` should be called exactly once at the end of this iteration round.