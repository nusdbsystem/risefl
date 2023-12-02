# Client side manual

need to `import di_zkp_interface`

The functionalities of a client is encapsulated in the python class `di_zkp_interface.ClientInterface`

A test run is in [test.py](../python/test.py).

## Creating a python class instance
To construct an instance of the client interface, run `di_zkp_interface.ClientInterface(num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
weight_bits, random_normal_bit_shifter,
num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
check_type, 
client_id, bul_pub_keys, bul_prv_key,
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
- `check_type`: the type of the check method. Supported:
  - TODO
  - TODO
- `client_id`: client ID. Note that client IDs start from 1. There is no client 0. 
- `bul_pub_keys`: vector of public keys of all the clients on the bulletin board. See [Public Bulletin Board](bulletin.md).
- `bul_prv_key`: this client's private key on the bulletin board.  See [Public Bulletin Board](bulletin.md).
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
`initialize_from_seed(seed)` initializes an FL training process. The string `seed` is used to generate independent group elements w<sub>j</sub>. It should be agreed among the server and all the clients. After initialization,we can start the iterative training process. For every iteration, we should call the following 5 methods:

## Class methods to call within every iteration (Suppose `client_id` is `i`)

### 1. Round 1: Generate a message that should be sent to the server
Given the l2 norm bound `norm_bound` and a vector of floats `u_float`, the class method `send_1(norm_bound, u_float)` returns a string that should be sent to the server in round 1. Internally, it also converts `u_float/norm_bound` (a vector of floats, all entries lie between -1 and 1) to fixed point integers. The return value encodes:
- the DHKE public key of client `i`, signed with `bul_prv_key`. (See [Public Bulletin Board](bulletin.md).)


`send_1(norm_bound, u_float)` should be called exactly once.

### 2. Round 2: Process a message received from the server, and generate a message that should be sent back to the server
`receive_and_send_2(bytes_str)` processes the string `bytes_str` which should be received from the server and returns a string that should be broadcast to all the clients in round 2. `bytes_str` should encode:
- the signed DHKE public keys from all the clients.

The class method checks the integrity of the signed DHKE public keys using `bul_pub_keys`, makes commitments to the weight updates and computes the verifiable Shamir's share of client `i`'s blind r<sub>i</sub>.

The return value encodes:
- the commitments of the weight updates of client `i`,
- the check string of client `i`'s blind r<sub>i</sub> and the DHKE-encrypted shares of r<sub>i</sub> via verifiable Shamir sharing.

`receive_and_send_2(bytes_str)` should be called exactly once.

### 3. Round 3: Process a message received from the server, and generate a message that should be sent back to the server
`receive_and_send_3(bytes_str)` processes the string `bytes_str` which should be received from the server and returns a string that should be broadcast to all the clients in round 3. `bytes_str` should encode:
- the check string of client `j`'s blind r<sub>j</sub> and  the DHKE-encrypted shares of client `j`'s blind r<sub>j</sub> for every `j != i`,
- the seed that generates the random normal samples,
- the corresponding multi-exponentiation values.

The class method checks the correctness of the multi-exponentiations, generates a proof based on the multi-exponentiations, decodes all the other clients' shares and checks their integrity.

The return value encodes:
- the proof that client `i`'s committed weight update has l2 norm at most `norm_bound` by probabilistic checking,
- the list of clients whose shares do not pass their check strings.

`receive_and_send_3(bytes_str)` should be called exactly once.

### 4. Round 4: Process a message received from the server, and generate a message that should be sent back to the server
`receive_and_send_4(bytes_str)` processes the string `bytes_str` which should be received from the server and returns a string that should be broadcast to all the clients in round 4. `bytes_str` should encode:
- the list of `j`s that report that client `i`'s share do not pass the check, provided that the list has length at most `max_malicious_clients`.

The return value encodes:
- the clear text shares of r<sub>i</sub> of the clients in the return value of `send_4(i)`, provided that the list has length at most `max_malicious_clients`.

`receive_and_send_4(bytes_str)` should be called exactly once.

### 5. Round 5: Process a message received from the server, and generate a message that should be sent back to the server
`receive_and_send_5(bytes_str)` processes the string `bytes_str` which should be received from the server and returns a string that should be broadcast to all the clients in round 5. `bytes_str` should encode:
- the clear text values of client `i`'s share of r<sub>j</sub> for `j`'s that don't pass the check strings,
- the list of malicious clients.

The return value encodes:
- the sum of client `i`'s shares of r<sub>j</sub> over all the valid client `j`s.

`receive_and_send_5(bytes_str)`  should be called exactly once.

