# from flask import request
# from flask import jsonify
# from flask import Flask
# from flask import Flask, render_template, request, jsonify
# app = Flask(__name__)
#
# # @app.route('/hello',methods=['POST'])
# # def hello():
# #     message = request.get_json(force=True)
# #     name = message['name']
# #     response = {
# #         'greeting': 'Hello, ' + name + '!'
# #     }
# #     return jsonify(response)
# @app.route('/rohan')
# def student():
#    # new_model = tf.keras.models.load_model("C:\\Users\\rohan\\PycharmProjects\\SciFaitEmory\\code\\saved_model.pb")
#    #
#    # # Check its architecture
#    # print(new_model.summary())
#
#    # return render_template('student.html')
#    return render_template("test.html",result = result)
#
# @app.route('/result',methods = ['POST', 'GET'])
# def result():
#    if request.method == 'POST':
#       result = request.form
#       return render_template("results.html",result = result)
#
# if __name__ == '__main__':
#    app.run(debug = True)
from flask import request
from flask import jsonify
from flask import Flask

app = Flask(__name__)

@app.route('/hello',methods=['POST'])
def hello():
    message = request.get_json(force=True)
    name = message['name']
    response = {
        'greeting': 'Hello, ' + name + '!'
    }
    return jsonify(response)

if __name__ == '__main__':
   app.run(debug = True)