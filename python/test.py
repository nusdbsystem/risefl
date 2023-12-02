import di_zkp_interface
import os
import base64

# initialize parameters

# number of clients
num_clients = 3

# maximum number of malicious clients
max_malicious_clients = 1

# number of parameters
dim = 2

# number of blinds per group element
num_blinds_per_group_element = 1

# number of bits of weight updates
weight_bits = 16

# random normal samples are multiplied by 2^random_normal_bit_shifter and rounded to the nearest integer
random_normal_bit_shifter = 24

# number of random normal samples
num_norm_bound_samples = 1000

# bit-bound of each inner product
inner_prod_bound_bits = weight_bits + random_normal_bit_shifter + 4

# bit-bound of sum of squares of inner products
max_bound_sq_bits = 2 * (weight_bits + random_normal_bit_shifter) + 20

# bound of l2-norm
norm_bound = 2.0

# a random string used to generate independent group elements, to be used by both the server and clients
random_bytes = os.urandom(64)
random_bytes_str = base64.b64encode(random_bytes).decode('ascii')
print("random_bytes_str = " + random_bytes_str)

# initialize server
server = di_zkp_interface.ServerInterface(
    num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
    weight_bits, random_normal_bit_shifter,
    num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
    di_zkp_interface.CHECK_TYPE_L2NORM,
    False)

server.initialize_from_seed(random_bytes_str)
server.initialize_new_iteration(norm_bound)

print("server.dim = " + str(server.dim))
print("server.weight_bits = " + str(server.weight_bits))

# create public keys on the bulletin board and the corresponding private keys
sign_pub_keys_vec = di_zkp_interface.VecSignPubKeys(num_clients + 1)
sign_prv_keys_vec = di_zkp_interface.VecSignPrvKeys(num_clients + 1)

for i in range(num_clients+1):
    sign_key_pair = di_zkp_interface.gen_sign_key_pair()
    sign_pub_keys_vec[i] = sign_key_pair.first
    sign_prv_keys_vec[i] = sign_key_pair.second

# initialize clients
clients = []
for i in range(num_clients+1):
    client = di_zkp_interface.ClientInterface(
        num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
        weight_bits, random_normal_bit_shifter,
        num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
        di_zkp_interface.CHECK_TYPE_L2NORM,
        i, sign_pub_keys_vec, sign_prv_keys_vec[i])
    # print("client.client_id = " + str(client.client_id))
    clients.append(client)

print("clients[0].dim = " + str(clients[0].dim))
print("clients[0].client_id = " + str(clients[0].client_id))

client_ids = range(1, num_clients+1)
print(client_ids)

for i in client_ids:
    clients[i].initialize_from_seed(random_bytes_str)

# iteration start:

# list of weight updates
# weight_updates_collection[0] is never used
# weight_updates_collection[i] is client i's weight updates, 1 <= i <= num_clients
weight_updates_collection = [[0, 0], [0, 1], [0.3, 0.4], [0.6, -0.9]]

# step 1: clients send messages to server
for i in client_ids:
    weights = di_zkp_interface.VecFloat(weight_updates_collection[i])
    server.receive_1(clients[i].send_1(norm_bound,
                                       weights), i)

print("***** step 1 finished *****")

# step 2: server sends (the same) message to clients, and gets messages back
bytes_sent_2 = server.send_2()
# print("bytes_sent_2 = " + bytes_sent_2)
for i in client_ids:
    client_i_str = clients[i].receive_and_send_2(bytes_sent_2)
    server.receive_2(client_i_str, i)

print("***** step 2 finished *****")

server.concurrent_process_before_send_3()

# step 3: server sends messages to clients, and gets messages back
for i in client_ids:
    server_send_3_str = server.send_3(i)
    client_i_str = clients[i].receive_and_send_3(server_send_3_str)
    server.receive_3(client_i_str, i)

print("***** step 3 finished *****")

# step 4: server sends messages to clients, and gets messages back
server.process_before_send_4()
for i in client_ids:
    server.receive_4(clients[i].receive_and_send_4(server.send_4(i)), i)

print("***** step 4 finished *****")

# step 5: server sends messages to clients, and gets messages back
server.process_before_send_5()
for i in client_ids:
    server.receive_5(clients[i].receive_and_send_5(server.send_5(i)), i)

print("***** step 5 finished *****")

# finish one iteration, aggregate sum is in server.final_update_float, average is in server.final_update_float_avg
server.finish_iteration()

col_sums = [sum(x) for x in zip(*weight_updates_collection)]
for j in range(dim):
    assert (abs(server.final_update_float[j] - col_sums[j]) < 1e-4)
    assert (abs(server.final_update_float[j] - num_clients * server.final_update_float_avg[j]) < 1e-4)

print("test python interface success")