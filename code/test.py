if __name__=="__main__":
    board = [["_" for i in range(3)] for j in range(3)]

    # define the players
    player1 = "X"
    player2 = "O"

    # define the current player
    current_player = player1

    # define the game state
    game_over = False

    while not game_over:
        # display the board
        for row in board:
            print(" ".join(row))

        # ask the current player for their move
        move = input(f"Player {current_player}, make your move (row, col): ")

        # parse the move
        row, col = map(int, move.split(","))

        # place the move on the board
        if board[row][col] == "_":
            board[row][col] = current_player
        else:
            print("Invalid move!")
            continue

        # check if the game is over
        if (
            # check rows
            all(board[row][0] == current_player for row in range(3))
            or all(board[row][1] == current_player for row in range(3))
            or all(board[row][2] == current_player for row in range(3))
            # check columns
            or all(board[0][col] == current_player for col in range(3))
            or all(board[1][col] == current_player for col in range(3))
            or all(board[2][col] == current_player for col in range(3))
            # check diagonals
            or all(board[i][i] == current_player for i in range(3))
            or all(board[i][2 - i] == current_player for i in range(3))
        ):
            print(f"Player {current_player} wins!")
            game_over = True
        elif all(board[row][col] != "_" for row in range(3) for col in range(3)):
            print("It's a tie!")
            game_over = True

        # switch to the other player
        if current_player == player1:
            current_player = player2
        else:
            current_player = player1
