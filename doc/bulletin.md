# Public Bulletin Board

The public bulletin board is used to ensure authenticity of the clients. Every client produces a key pair `(pub_key, prv_key)`. The public key `pub_key` is posted on the public bulletin board. The private key `prv_key` is used to sign messages by this client. Everyone can use `pub_key` to check authenticity of the signed messages. As long as a client trusts all the public keys of all the clients, the server can never steal information from the client by impersonation.

The function `di_zkp_interface.gen_sign_key_pair()` returns a pair of a public key (data type `di_zkp_interface.SignPubKey`) and the corresponding private key (data type `di_zkp_interface.SignPrvKey`). 