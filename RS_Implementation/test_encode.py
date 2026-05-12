from rs_codec.encoder import RSCodec, rs_encode_msg

def main():
    # Example message
    message = [32, 91, 11, 120, 209]
    nsym = 4  # number of parity symbols

    # Using helper function
    codeword = rs_encode_msg(message, nsym)
    print("\nEncoded codeword (from function):", codeword)

    # Using class
    codec = RSCodec(nsym)
    cw = codec.encode(message)
    print("Encoded codeword (from RSCodec):", cw)

if __name__ == "__main__":
    main()
