import os
from flask import Flask, flash, request, redirect, url_for, render_template
from werkzeug.utils import secure_filename, send_file
import tensorflow as tf
import predict

UPLOAD_FOLDER = '/Users/rohan/PycharmProjects/SciFaitEmory/code/uploads'
ALLOWED_EXTENSIONS = {'html','h5ad'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = b'_5#y2L"F4Q8z\n\xec]/'

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.errorhandler(413)
def too_large(e):
    return "File is too large", 413

@app.route("/cellcano")
def run_cellcano():
    return render_template("cellcano.html")

@app.route("/links")
def run_links():
    return render_template("links.html")
# @app.route('/download')
# def downloadFile ():
#     #For windows you need to use drive name [ex: F:/Example.pdf]
#     path = "/Examples.pdf"
#     return send_file(path, as_attachment=True)

@app.route('/upload', methods=['GET','POST'])
def upload_file():
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        #flash('file requested: '+file.filename)
        #flash('request URL'+request.url)
        # If the user does not select a file, the browser submits an
        # empty file without a filename.
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if not allowed_file(file.filename):
            flash('File type not supported, upload an h5ad file.')

        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            #flash('UPLOAD_FOLDER: '+app.config['UPLOAD_FOLDER'])
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

            # new_model = tf.keras.models.load_model('saved_model/my_model')
            # # now launch the NN from here, which should dump results in a file.
            # print(new_model)
            # # read results from the output file and set prediction_value to that...
            # prediction_value = type(new_model.summary())
            # print(prediction_value)


            #upload celltype file for download
            #completed above
                #use predict.py to generate file, search how to upload

            #return tsne plot
                # reference og code
            # marker gene heatmap too
                # ArchR or other R/python packages
                # target ot training?
            # data summary? proportions maybe
                # this can be done by parsing csv file
            #then delete used files once program executes
            test_adata = predict.predict()
            test_adata.obs[['pred_celltype']].to_csv(
                "C:\\Users\\rohan\\PycharmProjects\\SciFaitEmory\\code\\celltypes.csv")
            prediction_value = "celltypes.csv"
            return render_template("upload.html", prediction = prediction_value)

    return render_template("upload.html", prediction='No results/no upload')
if __name__ == "__main__":
    app.run(debug=True)